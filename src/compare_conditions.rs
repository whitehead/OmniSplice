#![allow(warnings)]

use std::collections::{HashMap, HashSet};
use std::fmt::{self, format};
use std::path::Path;
use std::sync::Arc;
use std::{fs, result};
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::time::Instant;

use clap::CommandFactory;
use clap::{Parser, Subcommand};
use nalgebra::coordinates::X;
use std::path::PathBuf;

use rayon::ThreadPoolBuilder;
use rayon::prelude::*;

//mod stat_common;
//use stat_common::errors::LogisticRegressionError;
//use stat_common::common::{Tester, parse_js_file, Genotype, JunctionStats, SplicingCategory,
//     welch_t_test, TtestResult, apply_bh_correction};
//use stat_common::glm_logistic::{GLM};
//use stat_common::glm_beta_binomial::GLMBetaBinomiale;

//use stat_common::common::{parse_js_file, Genotype, JunctionStats, SplicingCategory};
use std::process::{Command, Stdio};
mod common;

use common::error::OmniError;

use nalgebra::{DMatrix, DVector};
//use statrs::distribution::{ChiSquared, ContinuousCDF};
//use adjustp::{adjust, Procedure};

use crate::common::junction_file;
//use crate::stat_common::common::{TestResults, TestStatus};

use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};


use rayon::prelude::*;
use rayon::iter::Either;

///


 #[derive(Debug)]
pub struct CountsStats{
    spliced: u32,
    unspliced: u32,
    clipped: u32,
    exon_other: u32,
    skipped: u32,
    skipped_unrelated: u32,
    wrong_strand: u32,
    e_isoform: u32,
}

impl CountsStats{

    pub fn new(string: &[&str]) -> Self {
        let c = string
            .iter()
            .map(|x| x.parse::<u32>().unwrap())
            .collect::<Vec<u32>>();
        CountsStats {
            spliced: c[0],
            unspliced: c[1],
            clipped: c[2],
            exon_other: c[3],
            skipped: c[4],
            skipped_unrelated: c[5],
            wrong_strand: c[6],
            e_isoform: c[7],
        }
    }

    pub fn extract_(&self, category: &Vec<SplicingCategory>) -> u32{
        let mut count = 0;

        for cat in category{
            match cat{
                SplicingCategory::Spliced => {count += self.spliced}, 
                SplicingCategory::Unspliced => {count += self.unspliced},
                SplicingCategory::Clipped => {count += self.clipped},
                SplicingCategory::ExonOther => {count += self.exon_other},
                SplicingCategory::Skipped => {count += self.skipped}, 
                SplicingCategory::SkippedUnrelated => {count += self.skipped_unrelated}, 
                SplicingCategory::WrongStrand => {count += self.wrong_strand},
                SplicingCategory::EIsoform => {count += self.e_isoform}
            }
        }
        count
    }
}

#[derive(Debug)]
pub enum SplicingCategory{
    Spliced, 
    Unspliced,
    Clipped,
    ExonOther,
    Skipped, 
    SkippedUnrelated, 
    WrongStrand,
    EIsoform
}


impl TryFrom<&str> for SplicingCategory {
    type Error = &'static str;
    fn try_from(item: &str) -> Result<Self, Self::Error> {
        match item {
            "Spliced" | "SPLICED" => Ok(SplicingCategory::Spliced),
            "Unspliced" | "UNSPLICED" => Ok(SplicingCategory::Unspliced),
            "Clipped" | "CLIPPED" => Ok(SplicingCategory::Clipped),
            "Exon_other" | "EXONOTHER" => Ok(SplicingCategory::ExonOther),
            "Skipped" | "SKIPPED" => Ok(SplicingCategory::Skipped),
            "SkippedUnrelated"| "SKIPPEDUNRELATED" => Ok(SplicingCategory::SkippedUnrelated),
            "Wrong_strand" | "WRONGSTRAND" => Ok(SplicingCategory::WrongStrand),
            "E_isoform" | "EISOFORM" => Ok(SplicingCategory::EIsoform),
            _ => Err("input does not match one of the accepted values: Spliced Unspliced Clipped Exon_other Skipped SkippedUnrelated Wrong_strand  E_isoform")
        }
    }
}

impl fmt::Display for SplicingCategory {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        match self {
            SplicingCategory::Spliced => { write!(f, "SPLICED") }, 
            SplicingCategory::Unspliced => { write!(f, "UNSPLICED") },
            SplicingCategory::Clipped => { write!(f, "CLIPPED") },
            SplicingCategory::ExonOther => { write!(f, "EXONOTHER") },
            SplicingCategory::Skipped => { write!(f, "SKIPPED") }, 
            SplicingCategory::SkippedUnrelated => { write!(f, "SKIPPEDUNRELATED") }, 
            SplicingCategory::WrongStrand => { write!(f, "WRONGSTRAND") },
            SplicingCategory::EIsoform => { write!(f, "EISOFORM") }
        }
    }
}

#[derive(Debug)]
pub struct JunctionStats{
    pub contig: String, 
    pub start: String, 
    pub end: String, 
    pub strand: String,
    pub ambiguous: bool,
    pub control_count: Vec<CountsStats>,
    pub treat_count: Vec<CountsStats>,
    pub gene_tr: HashSet<String>,
    pub sample_done: HashSet<String>
}

impl JunctionStats{
    pub fn get_pos_string(&self) -> String{
        format!("{}:{}-{}({})", self.contig.clone(),
           self.start.clone(),
           self.end.clone(),
           self.strand.clone())
    }

   
}

impl Hash for JunctionStats{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.contig.hash(state); 
        self.strand.hash(state);
        self.start.hash(state); 
        self.end.hash(state);  
    }
}


fn header_to_map(header:  &str) -> Result<HashMap<String, usize>, OmniError>{
    Ok(header
        .trim()
        .split_whitespace()
        .enumerate()
        .map(|(i, v)| (v.to_string(), i))
        .collect())
    }

#[derive(Debug, PartialEq, Eq)]
pub enum Genotype{
    CONTROL,
    TREATMENT
}

pub fn parse_js_file(file_path: &str, result: &mut HashMap<String, JunctionStats>, genotype: Genotype) -> Result<(), OmniError>{



        let path = Path::new(file_path);
        let file_stem = format!("{:?}", path.file_stem().unwrap());
        info!("opening file");
        let f = File::open(file_path)?;
        let mut reader = BufReader::new(f);
        let mut header_line = String::new();

        info!("reading header");
        let len = reader.read_line(&mut header_line)?;
        let header = header_to_map(&header_line)?;

        
        let mut line: String = "".to_string();
        let mut spt: Vec<&str> = Vec::new();
        let mut trimed : &str = "";
        let mut gene: String = "".to_string();
        let mut contig: String = "".to_string();
        let mut strand : String = "".to_string();
        let mut start: String = "".to_string();
        let mut end: String = "".to_string();
        let mut key: String = "".to_string();

        for iterline in reader.lines(){
            line = iterline?;
            trimed = line.trim();

            if trimed.is_empty(){continue}
            spt = trimed.split("\t").collect::<Vec<&str>>();

            let ambi = if spt[*header.get("Ambiguous").unwrap()] == "true" {true} else {false}; 
  
            gene = format!("{}_{}_{}", spt[*header.get("Gene").unwrap()], spt[*header.get("Transcript").unwrap()], spt[*header.get("Intron").unwrap()]);
            contig = spt[*header.get("Contig").unwrap()].to_string();
            strand = spt[*header.get("Strand").unwrap()].to_string();
            start = if strand == "+" {spt[*header.get("Donnor").unwrap()].to_string()} else {spt[*header.get("Acceptor").unwrap()].to_string()};
            end = if strand == "+" {spt[*header.get("Acceptor").unwrap()].to_string()} else {spt[*header.get("Donnor").unwrap()].to_string()};
            
            key = format!("{} {} {} {}", contig, strand, start, end);


            // 3. case
            // never seen 
            // made in previous file
            // already seen in this file

            let mut never = false;
            let mut made = false; 
            let mut visited = false;

            if !result.contains_key(&key){
                never = true
            }
            else if result.get(&key).unwrap().sample_done.contains(&file_stem) {
                visited=true
            }
            else {
                made = true
            }



            if never{
                let counts = CountsStats::new(&spt[8..]);
                let mut junction = JunctionStats{
                        contig: contig, 
                        start: start, 
                        end: end, 
                        strand: strand,
                        ambiguous: ambi,
                        control_count: Vec::new(),
                        treat_count: Vec::new(),
                        gene_tr: HashSet::new(),
                        sample_done: HashSet::new()
                };
                junction.sample_done.insert(file_stem.clone());
                junction.gene_tr.insert(gene);
                match genotype {
                    Genotype::CONTROL => {junction.control_count.push(counts);}
                    Genotype::TREATMENT => {junction.treat_count.push(counts);}
                }
                result.insert(key.clone(), junction); // TODO remove this clonae after debug
            }
            else if made{
                let counts = CountsStats::new(&spt[8..]);
                match result.get_mut(&key){
                    Some(j) => {
                        j.sample_done.insert(file_stem.clone());
                        j.gene_tr.insert(gene);
                        match genotype {
                        Genotype::CONTROL => {j.control_count.push(counts);}
                        Genotype::TREATMENT => {j.treat_count.push(counts);}
                        }
                    },
                     None => {warn!("unreachable 'made' case in parse js"); ()}
                }
                
            }
            else if visited{
                match result.get_mut(&key){
                    Some(j) => {
                        j.sample_done.insert(file_stem.clone());
                        j.gene_tr.insert(gene);
                    },
                    None => {warn!("unreachable 'visited' case in parse js"); ()}
            }
        }

        /*if spt[*header.get("Gene").unwrap()] == "FBgn0000015"{
                warn!("{:?}", spt);
                warn!("{:?} {:?} {:?} {:?}", key.clone(), never, visited, made);
                warn!("{:?}", result.get(&key));

            } */

    }
    Ok(())
}



///




fn run_R(in_file: &str, out_file: &str, thread: u32,
     min_read:u32, min_unsp:u32, do_ambi: bool, no_beta: bool){
    let script_path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("R/R_osStat.R");
    
    let mut bowt_child = Command::new("Rscript")
        .args([script_path.to_str().unwrap(), in_file, out_file,
          &thread.to_string(), &min_read.to_string(),
          &min_unsp.to_string(), &no_beta.to_string(), &do_ambi.to_string()])
        .spawn()
        .expect("Rscript failed to start");
    let _result = bowt_child.wait().unwrap();
}

//}
// If this ahppend to be to expansive I could order the indice of qvalue and then retorder both list accordingly.
fn sort_by_f32_copy<T>(data: &mut Vec<T>, scores: &mut Vec<f32>) {
    // Create pairs, sort them, then unzip
    let mut pairs: Vec<(T, f32)> = data.drain(..).zip(scores.drain(..)).collect();
    
    pairs.sort_by(|a, b| {
        a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal)
    });
    
    let (sorted_data, sorted_scores): (Vec<T>, Vec<f32>) = pairs.into_iter().unzip();
    
    *data = sorted_data;
    *scores = sorted_scores;
}



fn sort_by_f64_permutation<T: std::fmt::Debug>(data: &mut Vec<T>,
                                               scores: &mut Vec<f64>) {

    let mut indices: Vec<usize> =  (0..scores.len()).collect();
    let mut rank = vec![0usize; scores.len()];
    {
        indices.sort_by(|&a, &b| scores[a]
            .partial_cmp(&scores[b])
            .unwrap_or_else(|| std::cmp::Ordering::Equal));
        
        
        for (sorted_pos, &orig_idx) in indices.iter().enumerate() {
            rank[orig_idx] = sorted_pos;
        }
    }
        
    let n = rank.len();
    
        for i in 0..n{
        let mut j = rank[i];
        while j != i{
            data.swap(i, j);
            scores.swap(i, j);
            rank.swap(i, j);
            j = rank[i];
        }
    }

}





fn helper_6(value: Option<f64>) -> String{
    match value{
        Some(v) => format!("{:e}", v),
        None => "nan".to_string()
    }
}

//["chr", "strand", "start", "end",
 //
 //   "control_success", "control_failures", "control_ratio",
 //   "treatment_success", "treatment_failures", "treatment_ratio",
 //   "glmbb_p_value", "glmbb_q_value", "oddRatio", "glmbb_CI", "glmbb_status", 
 //   "ttest_stat", "ttest_pvalue", "ttest_q_value", 
 //   "gene_transcript_intron"];




fn get_cohensh(pcontrol: f32, ptreat: f32) -> f32{
    2.0 * (ptreat.sqrt().asin() - pcontrol.sqrt().asin()) //f64.asin() 
}
/* 
fn parse_results_update_vec(vec_r: &Vec<(&JunctionStats, TestResults, TtestResult)>,
                            result: &mut Vec<Vec<String>>){

    let mut value: (String, String, String, String);
    for (i, (j, t, tt)) in vec_r.into_iter().enumerate(){
        
        let mut f: Vec<String> = Vec::new();
        f.push(j.get_pos_string());

        value = t.string_count.clone();

        f.push(value.0);
        f.push(value.1);
        f.push(   match t.control_prop {
                        Some(c) => c.to_string(),
                        None => "nan".to_string()}
        );

        f.push(value.2);
        f.push(value.3);
        f.push(  match t.treatment_prop {
            Some(c) => c.to_string(),
            None => "nan".to_string(),}
        );

        match (t.control_prop, t.treatment_prop){
            (Some(p_c), Some(p_t)) => {
                f.push(get_cohensh(p_c, p_t).to_string());
                f.push((p_t - p_c).to_string());
            }
            _ => {
                f.push("nan".to_string());
                f.push("nan".to_string());
            }
        }
//"glmbb_p_value", "glmbb_q_value", "oddRatio", "glmbb_CI", "glmbb_status", 

        f.push(helper_6(t.p_value));
        f.push(helper_6(t.q_value));

        match  t.odd_ratio{ 
            Some(odd) => f.push(format!("{:6.2e}", odd)),
            None => f.push("nan".to_string())
        }

        match  (t.or_ci_lower, t.or_ci_upper){ 
            (Some(low), Some(up)) => f.push(format!("{:4.2e}/{:4.2e}",low, up)),
            _ => f.push("nan/nan".to_string())
        }

        match &t.status{
            Some(status) =>  f.push(format!("{}", status)),
            None => f.push("nan".to_string())
        }


        match tt.t_stat {
            Some(r) => f.push(format!("{:6.2e}", r)),
            None => f.push("nan".to_string())

        }

        f.push(helper_6(tt.p_value));
        f.push(helper_6(tt.q_value));



        f.push(j.gene_tr.iter().map(|x| x.to_owned()).collect::<Vec<String>>().join(";"));
        result.push(f)
    }
}
*/



fn get_count(junction: &JunctionStats, 
                        successes_cat: &Vec<SplicingCategory>,           
                        failures_cat: &Vec<SplicingCategory>) -> (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>){

        let mut ctrl_suc: Vec<u32> = Vec::new();
        let mut ctrl_fail: Vec<u32> = Vec::new();
        let mut treat_suc: Vec<u32> = Vec::new();
        let mut treat_fail: Vec<u32> = Vec::new();

        for count in &junction.control_count{
            ctrl_suc.push(count.extract_(&successes_cat));
            ctrl_fail.push(count.extract_(&failures_cat));     
        }
        for count in &junction.treat_count{
            treat_suc.push(count.extract_(&successes_cat));
            treat_fail.push(count.extract_(&failures_cat));     
        }
        
        (ctrl_suc, ctrl_fail, treat_suc, treat_fail)
    }


pub fn pass_min_read(min_cover: u32, min_unspliced: u32, counts: (&Vec<u32>, &Vec<u32>, &Vec<u32>, &Vec<u32>)) -> bool {
        let (g1_succ, g1_fail, g2_succ, g2_fail) = counts;
        if (g1_fail.iter().sum::<u32>() + g1_succ.iter().sum::<u32>() < min_cover) ||
                 (g2_fail.iter().sum::<u32>() + g2_succ.iter().sum::<u32>() < min_cover){
            return false;
        }
        if g1_fail.iter().sum::<u32>() < min_unspliced && g2_fail.iter().sum::<u32>() < min_unspliced{
            return false;
        }

        if g1_succ.iter().sum::<u32>() < min_unspliced && g2_succ.iter().sum::<u32>() < min_unspliced{
            return false;
        }
        true
}

fn run_one_test(junction: &HashMap<String, JunctionStats>,
                successes_cat: Vec<SplicingCategory>,
                failures_cat: Vec<SplicingCategory>,
                out_file_path: &str,
                ambi: bool, min_cover: u32, min_fail:u32,
                thread:u32, low_rep: bool, control: Vec<String>, treatment: Vec<String>) -> Result<(), Box<dyn std::error::Error + Send + Sync>>{
     let successes_cat_ref = Arc::new(successes_cat);
     let failures_cat_ref = Arc::new(failures_cat);
    
    println!("Starting test!");
    println!("Starting compiling count!");
        let now = Instant::now();



    let vec_res = junction
        .par_iter()
        .filter_map(|(k, j)| {

            let (ctrl_suc, ctrl_fail, treat_suc, treat_fail) =
             get_count(&j, &successes_cat_ref, &failures_cat_ref);


            if !pass_min_read(min_cover, min_fail,(&ctrl_suc, &ctrl_fail, &treat_suc, &treat_fail)){return  None}
            if j.ambiguous == true && ambi == true {return None}

            let mut f= vec![j.get_pos_string()];
            f.push(if j.ambiguous {"true".to_string()} else {"false".to_string()});


            let ctrl_suc_sum = ctrl_suc.iter().fold(0, |acc, x| acc + x);
            let ctrl_fail_sum = ctrl_fail.iter().fold(0, |acc, x| acc + x);
            let ctrl_prop = if (ctrl_suc_sum + ctrl_fail_sum) > 0 { Some(ctrl_suc_sum as f32 / (ctrl_suc_sum as f32  + ctrl_fail_sum as f32 ))} else {None};

            f.push(ctrl_suc.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
            f.push(ctrl_fail.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));

            f.push(match ctrl_prop  {
                          Some(c) => c.to_string(),
                         None => "nan".to_string()}
            );

            let treat_suc_sum = treat_suc.iter().fold(0, |acc, x| acc + x);
            let treat_fail_sum = treat_fail.iter().fold(0, |acc, x| acc + x);
            let treat_prop = if (treat_suc_sum + treat_fail_sum) > 0 { Some(treat_suc_sum as f32  / (treat_fail_sum as f32  + treat_suc_sum as f32 ))} else {None};
  
            f.push(treat_suc.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
            f.push(treat_fail.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","));
            f.push( match treat_prop{
                Some(c) => c.to_string(),
                None => "nan".to_string(),}
            );


            match (ctrl_prop, treat_prop){
                (Some(p_c), Some(p_t)) => {
                    f.push(get_cohensh(1.0 - p_c, 1.0 - p_t).to_string());
                    f.push(((1.0-p_t) - (1.0 - p_c)).to_string());
                }
                _ => {
                    f.push("nan".to_string());
                    f.push("nan".to_string());
                }
            }
            f.push(j.gene_tr.iter().map(|x| x.to_owned()).collect::<Vec<String>>().join(";"));
        Some(f.join("\t"))
    }).collect::<Vec<String>>();

    let elapsed_time = now.elapsed();
    println!("compiling counts tooks {} seconds.", elapsed_time.as_secs());   

    let mut out_file_open =
        File::create_new(out_file_path.clone()) //presorted out_file.clone()
            .unwrap_or_else(|_| panic!("output file {} should not exist.", &out_file_path)); //expect(&format!("output file {} should not exist.", &table));
    let mut out_stream = BufWriter::new(out_file_open);


    out_stream.write(format!("#CONTROL_FILE: {}; TREATMENT_FILE: {}\n", control.join(","), treatment.join(",")).as_bytes());
    out_stream.write(format!("#success: {} ; failures: {}\n",
        successes_cat_ref.iter().map(|x| format!("{}", x).to_string()).collect::<Vec<String>>().join(" "),
        failures_cat_ref.iter().map(|x| format!("{}", x).to_string()).collect::<Vec<String>>().join(" ").to_string()).as_bytes());
    out_stream.write(format!("#min_read: {}; min_fail {}; ambigious: {}, test: {}\n", min_cover, min_fail, ambi, if low_rep == true{"logistic/Fischer"}else{"beta-binomial/t-test"}).as_bytes());
    let header = vec!["chr:start-end(strand)", "ambigious",
    "control_success", "control_failures", "control_ratio",
     "treatment_success", "treatment_failures", "treatment_ratio",  "Cohens-h", "delta-psi",
       "gene_transcript_intron"];

    out_stream.write(format!("{}\n", header.join("\t")).as_bytes());
    for e in vec_res{
        out_stream.write(format!("{}\n", e).as_bytes())?;
    }

    let _ = out_stream.flush();


    // TODO add ambigious 
    // TODO add do beta option.

    let elapsed_time = now.elapsed();
let now = Instant::now();   
    run_R(out_file_path, out_file_path, thread,
         min_cover, min_fail,
        ambi, low_rep);
    println!("Running R tooks {} seconds.", elapsed_time.as_secs());   
    

    Ok(())
}




#[derive(Parser)]
#[command(name = "compare")]
#[command(about = "Allows to compare condition from Omnisplice junction file", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Run a single comparison
    Run {
        /// Output file prefix path
        #[arg(short, long,  required = true)]
        outfile: PathBuf,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        control_files: Vec<PathBuf>,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        treatment_files: Vec<PathBuf>,

        /// Splicing type considered as good
        /// must be one or more of the follwing: Spliced Unspliced  Clipped Exon_other Skipped SkippedUnrelated Wrong_strand  E_isoform
        #[arg(short, long, num_args = 1.., required = true)]
        splicing_ok: Vec<String>,

        /// Splicing type considered as failures
        /// must be one or more of the follwing: Spliced Unspliced       Clipped Exon_other      Skipped SkippedUnrelated        Wrong_strand    E_isoform
        #[arg(short, long, num_args = 1.., required = true)]
        splicing_fail: Vec<String>,

       /// Do you want to consider ambigious junction (overlaping exon)
       #[arg( long,)]
       ambigious: bool,

       /// Minimum read count for test succes + failure in each group must be >= min_read. Discards junctions (p-value = NaN) .
       /// Default: 30 |
       #[arg( long, default_value_t = 30)]
       min_read: u32,


       /// The data is organize as the count for 4 categories: control succ, control fail, treat succ, treat fail
       /// if  control succ < min_fail AND control fail < min_fail  -> fail the the test
       /// if  treat succ < min_fail AND treat fail < min_fail  -> fail the the test
       /// Default: 10
       #[arg( long, default_value_t = 10)]
       min_fail: u32,

    /// thread number default 5
    #[arg( long, default_value_t = 5)]
    thread: usize,

    /// by default OmniSplice uses glm logit / fischer test.
    /// toggle this flag to use a beta binomial and t-test instead. 
    #[arg(long, action = clap::ArgAction::SetFalse)]
    beta_bin: bool


    },
    /// Run all comparisons against splices
    RunAll {
        /// Output file path
        #[arg(short, long, required = true)]
        outfile_prefix: String,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        control_files: Vec<PathBuf>,

        /// Control Condition
        #[arg(short, long, num_args = 1.., required = true)]
        treatment_files: Vec<PathBuf>,

        /// The data is organize as the count for 4 categories: control succ, control fail, treat succ, treat fail
       /// if  control succ < min_read OR treat succ < min_read  -> fail the the test
       /// Default: 30
       #[arg( long, default_value_t = 30)]
       min_read: u32,

       /// The data is organize as the count for 4 categories: control succ, control fail, treat succ, treat fail
       /// if  control succ < min_fail AND control fail < min_fail  -> fail the the test
       /// if  treat succ < min_fail AND treat fail < min_fail  -> fail the the test
       /// Default: 10
       #[arg( long, default_value_t = 10)]
       min_fail: u32,

              /// thread number default 5
        #[arg( long, default_value_t = 5)]
        thread: usize,

    /// by default OmniSplice uses glm logit / fischer test.
    /// toggle this flag to use a beta binomial and t-test instead. 
    #[arg(long, action = clap::ArgAction::SetFalse)]
    beta_bin: bool
    },
}

fn parse_cat(input: Vec<String>) -> Result<Vec<SplicingCategory>, &'static str >{
    let mut res = Vec::new();
    for e in input{
        if e.trim().is_empty(){
            continue
        }
        res.push(SplicingCategory::try_from(e.trim())?)
    }
    Ok(res)
}

fn main() -> Result<(), Box<dyn std::error::Error + Send + Sync>> {


    Logger::try_with_str("error").unwrap().start().unwrap();


    let cli = Cli::parse();
    //env_logger::init();

    match cli.command {
        Commands::Run { outfile, control_files,
                        treatment_files, splicing_ok,
                        splicing_fail, ambigious, min_read , min_fail, 
                    thread, beta_bin} => {
        
            info!("Running comparison, output: {:?}", outfile);

                assert!(thread > 0);
            rayon::ThreadPoolBuilder::new()
            .num_threads(thread)          // your cap
            .build_global()
            .unwrap();

            let control = parse_cat(splicing_ok)?;
            let treatment = parse_cat(splicing_fail)?;
    
            let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);

            for file in &control_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::CONTROL).unwrap();
                    info!("done reading");
                }

            for file in &treatment_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::TREATMENT).unwrap();
                    info!("done reading");
                }


        run_one_test( &res, control,
                 treatment,  outfile.to_str().unwrap(), 
                  ambigious, min_read,
                   min_fail, thread as u32, beta_bin,
                    control_files.into_iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.into_iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;

        }
        Commands::RunAll { control_files, 
            treatment_files, outfile_prefix, 
            min_read , min_fail, thread, beta_bin} => {
            println!("Running all single comparisons, output: {:?}", outfile_prefix);
            println!("using logistic model and Fisher test: {}", !beta_bin);
            assert!(thread > 0, "thread must be postive");

            rayon::ThreadPoolBuilder::new()
            .num_threads(thread)          
            .build_global()
            .unwrap();

            let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);

        let now = Instant::now();

            for file in &control_files.clone(){
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::CONTROL).unwrap();
                    info!("done reading");
                }

            for file in &treatment_files.clone(){
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::TREATMENT).unwrap();
                    info!("done reading");
                }
            let elapsed_time = now.elapsed();
             println!("reading junction files tooks {} seconds.", elapsed_time.as_secs());   

    println!("All junction file parsed");

    let shared = Arc::new(res);

    let mut p = Path::new(&outfile_prefix).to_path_buf();


    //let mut jobs: Vec<(Vec<SplicingCategory>, Vec<SplicingCategory>, &str, bool, u32, u32)> = Vec::new();
    let mut p = Path::new(&outfile_prefix).to_path_buf();

    info!( "starting Unspliced" );
    let _ = p.set_extension("Unspliced.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                 vec![SplicingCategory::Unspliced],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, beta_bin,
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;



    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Unspliced],
    //   p.to_str().unwrap() , false, min_read, min_fail));
    info!( "starting WrongStrand" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("WrongStrand.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::WrongStrand],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, beta_bin, 
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;

    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::WrongStrand],
    //   p.to_str().unwrap() , false, min_read, min_fail));
    info!( "starting Skipped" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Skipped.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                 vec![SplicingCategory::Skipped],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, beta_bin,
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Skipped],
    //   p.to_str().unwrap() , false, min_read, min_fail));
        info!( "starting SkippedUnrelated" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("SkippedUnrelated.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::SkippedUnrelated],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, beta_bin,
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::SkippedUnrelated],
    //   p.to_str().unwrap() , false, min_read, min_fail));
        info!( "starting Clipped" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Clipped.tsv");

    run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::Clipped],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, beta_bin,
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;
   
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Clipped],
    //   p.to_str().unwrap() , false, min_read, min_fail));
    info!( "starting Exon_other" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("ExonOther.tsv");

        run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::ExonOther],  p.to_str().unwrap(), 
                  true, min_read, min_fail, thread as u32, beta_bin,
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::ExonOther],
    //   p.to_str().unwrap(), true, min_read, min_fail));
           info!( "starting E_isoform" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Isoform.tsv");
            run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::EIsoform],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, beta_bin,
                control_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>(),
                     treatment_files.iter().map(|x| x.display().to_string()).collect::<Vec<String>>())?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::EIsoform],
    //   p.to_str().unwrap() , false, min_read, min_fail));


    //println!("{:?}", );
    let now = Instant::now();
    let elapsed_time = now.elapsed();
    println!("Running slow_function() took {} seconds.", elapsed_time.as_secs());   

        }
    }

    Ok(())
}



       
