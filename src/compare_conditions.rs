#![allow(warnings)]

use std::collections::{HashMap, HashSet};
use std::fmt::format;
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

mod stat_common;
use stat_common::errors::LogisticRegressionError;
use stat_common::common::{Tester, parse_js_file, Genotype, JunctionStats, SplicingCategory,
     welch_t_test, TtestResult, apply_bh_correction};
use stat_common::glm_logistic::{GLM};
use stat_common::glm_beta_binomial::GLMBetaBinomiale;

use std::process::{Command, Stdio};
mod common;

use common::error::OmniError;

use nalgebra::{DMatrix, DVector};
use statrs::distribution::{ChiSquared, ContinuousCDF};
use adjustp::{adjust, Procedure};

use crate::common::junction_file;
use crate::stat_common::common::{TestResults, TestStatus};

use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};


use rayon::prelude::*;
use rayon::iter::Either;



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


fn run_one_test(junction: &HashMap<String, JunctionStats>,
                successes_cat: Vec<SplicingCategory>,
                failures_cat: Vec<SplicingCategory>,
                out_file_path: &str,
                ambi: bool, min_cover: u32, min_fail:u32,
                thread:u32, no_beta: bool) -> Result<(), Box<dyn std::error::Error + Send + Sync>>{

    
    info!("Starting test!");
    let vec_res = junction
        .par_iter()
        .map(|(k, j)| {
            let mut glm = GLMBetaBinomiale::new(
                &j.treat_count,
                &j.control_count,
                &successes_cat,
                &failures_cat,
                k.to_owned());
            let mut f= vec![j.get_pos_string()];
            f.push(if j.ambiguous {"true".to_string()} else {"false".to_string()});

            let t = glm.test(false, min_cover, min_fail);


            let value = t.string_count.clone();

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
            f.push(j.gene_tr.iter().map(|x| x.to_owned()).collect::<Vec<String>>().join(";"));
        f.join("\t")
    }).collect::<Vec<String>>();

    
    /*let (mut vec_ok, mut vec_err): (Vec<_>, Vec<_>) = junction

        .par_iter()
        .map(|(k, j)| {
            let mut glm = GLMBetaBinomiale::new(
                &j.treat_count,
                &j.control_count,
                &successes_cat,
                &failures_cat,
                k.to_owned(),
            );

            let x = if ambi == false && j.ambiguous == true{
                glm.test(true, min_cover, min_fail)
            }
            else {
                glm.test(false, min_cover, min_fail)
            };

            match x.status {
                Some(TestStatus::TreatmentIsNull)|
                Some(TestStatus::ControlIsNull) | Some(TestStatus::FailFilter) => {
                        Either::Right((j, x, TtestResult::new_empty()))
                        //vec_err.push((&j, x, TtestResult::new_empty()));
                },

                _ => {
                    let (g1_succ, g1_fail, g2_succ, g2_fail) = glm.get_sub_vec();
                    let g1_prop = g1_succ.iter().zip(g1_fail.iter()).map(|(x, y)| *x as f64  / (*x as f64 + *y as f64)).collect::<Vec<f64>>();
                    let g2_prop = g2_succ.iter().zip(g2_fail.iter()).map(|(x, y)| *x as f64  / (*x as f64 + *y as f64)).collect::<Vec<f64>>();

                    let t_testR = welch_t_test(&g1_prop, &g2_prop).unwrap_or_else(|_| TtestResult::new_empty());
                    
                    Either::Left((j, x, t_testR))
                }
            }
            
        }).partition_map(|either| either);
    /*let mut vec_ok: Vec<(&JunctionStats, TestResults, TtestResult)> = Vec::new();
    let mut vec_err: Vec<(&JunctionStats, TestResults, TtestResult)> = Vec::new();
    let mut x: TestResults; 
    for (k, j ) in junction.iter(){
        let mut glm = GLMBetaBinomiale::new(  
                                        &j.treat_count,
                                        &j.control_count,                                      
                                                &successes_cat,
                                                &failures_cat,
                                            k.to_owned());


        if ambi == false && j.ambiguous == true{
            x = glm.test(true, min_cover, min_fail);
        }
        else {
            x = glm.test(false, min_cover, min_fail);
        }

         match x.status {
            Some(TestStatus::TreatmentIsNull) | Some(TestStatus::ControlIsNull) |
            Some(TestStatus::ControlIsNull) | Some(TestStatus::FailFilter) => {

                    vec_err.push((&j, x, TtestResult::new_empty()));
            },

            _ => {
                let (g1_succ, g1_fail, g2_succ, g2_fail) = glm.get_sub_vec();
                let g1_prop = g1_succ.iter().zip(g1_fail.iter()).map(|(x, y)| *x as f64  / (*x as f64 + *y as f64)).collect::<Vec<f64>>();
                let g2_prop = g2_succ.iter().zip(g2_fail.iter()).map(|(x, y)| *x as f64  / (*x as f64 + *y as f64)).collect::<Vec<f64>>();

                let t_testR = match welch_t_test(&g1_prop, &g2_prop){
                    Ok(t) => t,
                    _ => TtestResult::new_empty()
                }; 
                
                vec_ok.push((&j, x, t_testR));
            }
         }

        /*if x.p_value.is_some(){

        }

        match x.status {
            Some(TestStatus::TreatmentIsNull) | Some(TestStatus::ControlIsNull) |
            Some(TestStatus::ControlIsNull) | Some(TestStatus::FailFilter) => {
                    vec_err.push((&j, x));
            },


        }

        let (g1_succ, g1_fail, g2_succ, g2_fail) = glm.get_sub_vec();
        let g1_prop = g1_succ.iter().zip(g1_fail.iter()).map(|(x, y)| *x as f64  / (*x as f64 + *y as f64)).collect::<Vec<f64>>();
        let g2_prop = g2_succ.iter().zip(g2_fail.iter()).map(|(x, y)| *x as f64  / (*x as f64 + *y as f64)).collect::<Vec<f64>>();

        let t_testR = welch_t_test(&g1_prop, &g2_prop);
        
        if x.p_value.is_some(){
            vec_ok.push((&j, x));
        }
        else{
            match x.status {
                Some(TestStatus::TreatmentIsNull) | Some(TestStatus::ControlIsNull) | Some(TestStatus::ControlIsNull) | Some(TestStatus::FailFilter) => {
                    vec_err.push((&j, x));
                } ,
                _ => vec_err.push((&j, x)),
                None => ()
            }

        }*/

    }*/

    apply_bh_correction(&mut vec_ok,
    |x| x.1.p_value.map(|p| p as f64),
    |x, q| x.1.q_value = Some(q));

    apply_bh_correction(&mut vec_ok,
    |x| x.2.p_value.map(|p| p as f64),
    |x, q| x.2.q_value = Some(q));


    /*let (mut vec_ok, mut vec_err): 
        (Vec<(&JunctionStats, TestResults, TtestResult)>,
         Vec<(&JunctionStats, TestResults, TtestResult)>) = vec_ok.into_iter().partition(|x| x.1.p_value.is_some());

    let pval = vec_ok.iter().map(|x| x.1.p_value.unwrap() as f64).collect::<Vec<f64>>();
    let mut qvalues = adjust(&pval, Procedure::BenjaminiHochberg);
    vec_ok.iter_mut().zip(qvalues.into_iter()).for_each(|(elem, qval)| {
        elem.1.q_value = Some(qval);
    });

    vec_ok.append(&mut vec_err);

        let (mut vec_ok, mut vec_err): 
        (Vec<(&JunctionStats, TestResults, TtestResult)>,
         Vec<(&JunctionStats, TestResults, TtestResult)>) = vec_ok.into_iter().partition(|x| x.1.p_value.is_some());
    */


    
    let mut final_: Vec<Vec<String>> = Vec::new();
    let mut value: (String, String, String, String);

    vec_ok.sort_by(|a, b|
        a.1.q_value.unwrap_or(10.).partial_cmp(
            &b.1.q_value.unwrap_or(10.)).unwrap()
        );


    parse_results_update_vec(&vec_ok, &mut final_);
    //parse_results_update_vec(&vec_err, &mut final_, None);

*/
    
    //let mut out_file_path_temp = Path::new(&out_file_path).to_path_buf();
   // let _ = out_file_path_temp.set_extension(".tmp");

    let mut out_file_open =
        File::create_new(out_file_path.clone()) //presorted out_file.clone()
            .unwrap_or_else(|_| panic!("output file {} should not exist.", &out_file_path)); //expect(&format!("output file {} should not exist.", &table));
    let mut out_stream = BufWriter::new(out_file_open);


    out_stream.write(format!("#success: {}\n", successes_cat.into_iter().map(|x| format!("{}", x).to_string()).collect::<Vec<String>>().join(" ")).to_string().as_bytes());
    out_stream.write(format!("#failures: {}\n", failures_cat.into_iter().map(|x| format!("{}", x).to_string()).collect::<Vec<String>>().join(" ")).to_string().as_bytes());
    out_stream.write(format!("#min_read: {}; min_fail {}; ambigious: {}\n", min_cover, min_fail, ambi).to_string().as_bytes());
    let header = vec!["chr:start-end(strand)", "ambigious",
    "control_success", "control_failures", "control_ratio",
     "treatment_success", "treatment_failures", "treatment_ratio", "delta-psi", "Cohens-h",
       "gene_transcript_intron"];//"glmbb_p_value", "glmbb_q_value", "oddRatio", "glmbb_CI", "glmbb_status", 
        //"ttest_stat", "ttest_pvalue", "ttest_q_value", 
        //   "gene_transcript_intron"];

    out_stream.write(format!("{}\n", header.join("\t")).as_bytes());
    for e in vec_res{
        out_stream.write(format!("{}\n", e).as_bytes())?;
    }

    let _ = out_stream.flush();


    // TODO add ambigious 
    // TODO add do beta option.
    run_R(out_file_path, out_file_path, thread,
         min_cover, min_fail,
        ambi, no_beta);

    

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

        /// Minimum unspliced read count for test  failure in at least one group must be >= min_fail. Discards junctions (p-value = NaN) .
       /// Default: 10 | 
       #[arg( long, default_value_t = 10)]
       min_fail: u32,

    // thread number default 5
    #[arg( long, default_value_t = 5)]
    thread: usize,

    // deactivate glm_bb (usefull if time is an issue)
    #[arg(long, action = clap::ArgAction::SetFalse)]
    no_beta: bool


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

       /// Minimum read count for test succes + failure in each group must be >= min_read. Discards junctions (p-value = NaN) .
       /// Default: 30
       #[arg( long, default_value_t = 30)]
       min_read: u32,

        /// Minimum unspliced read count for test  failure in at least one group must be >= min_fail. Discards junctions (p-value = NaN) .
       /// Default: 10
       #[arg( long, default_value_t = 10)]
       min_fail: u32,

              // thread number default 5
        #[arg( long, default_value_t = 5)]
        thread: usize,

        // deactivate glm_bb (usefull if time is an issue)
        #[arg(long, action = clap::ArgAction::SetFalse)]
        no_beta: bool
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
                    thread, no_beta} => {
        
            info!("Running comparison, output: {:?}", outfile);

                assert!(thread > 0);
            rayon::ThreadPoolBuilder::new()
            .num_threads(thread)          // your cap
            .build_global()
            .unwrap();

            let control = parse_cat(splicing_ok)?;
            let treatment = parse_cat(splicing_fail)?;
    
            let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);

            for file in control_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::CONTROL).unwrap();
                    info!("done reading");
                }

            for file in treatment_files{
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::TREATMENT).unwrap();
                    info!("done reading");
                }


        run_one_test( &res, control,
                 treatment,  outfile.to_str().unwrap(), 
                  ambigious, min_read,
                   min_fail, thread as u32, no_beta)?;

        }
        Commands::RunAll { control_files, 
            treatment_files, outfile_prefix, 
            min_read , min_fail, thread, no_beta} => {
            info!("Running all single comparisons, output: {:?}", outfile_prefix);
                
            assert!(thread > 0, "thread must be postive");
            rayon::ThreadPoolBuilder::new()
            .num_threads(thread)          // your cap
            .build_global()
            .unwrap();
            let mut res:  HashMap<String, JunctionStats> = HashMap::with_capacity(1_000_000);

            for file in control_files.clone(){
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::CONTROL).unwrap();
                    info!("done reading");
                }

            for file in treatment_files.clone(){
                    info!("parsing {:?}", file);
                    parse_js_file(file.to_str().unwrap(), &mut res, Genotype::TREATMENT).unwrap();
                    info!("done reading");
                }
            
            //println!("{:?}", res.get("3R - 16931613 16935926"));
            //println!("{:?}", &treatment_files);
            //println!("{:?}", &control_files);

       //     ThreadPoolBuilder::new()
       // .num_threads(6)          // ← set the limit here
       // .build_global()
       // .expect("Failed to initialise Rayon thread‑pool");

    info!("All junction file parsed");

    let shared = Arc::new(res);

    let mut p = Path::new(&outfile_prefix).to_path_buf();


    //let mut jobs: Vec<(Vec<SplicingCategory>, Vec<SplicingCategory>, &str, bool, u32, u32)> = Vec::new();
    let mut p = Path::new(&outfile_prefix).to_path_buf();

    info!( "starting Unspliced" );
    let _ = p.set_extension("Unspliced.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                 vec![SplicingCategory::Unspliced],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, no_beta)?;



    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Unspliced],
    //   p.to_str().unwrap() , false, min_read, min_fail));
    info!( "starting WrongStrand" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("WrongStrand.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::WrongStrand],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, no_beta)?;

    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::WrongStrand],
    //   p.to_str().unwrap() , false, min_read, min_fail));
    info!( "starting Skipped" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Skipped.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                 vec![SplicingCategory::Skipped],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, no_beta)?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Skipped],
    //   p.to_str().unwrap() , false, min_read, min_fail));
        info!( "starting SkippedUnrelated" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("SkippedUnrelated.tsv");
    run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::SkippedUnrelated],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, no_beta)?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::SkippedUnrelated],
    //   p.to_str().unwrap() , false, min_read, min_fail));
        info!( "starting Clipped" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Clipped.tsv");

    run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::Clipped],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, no_beta)?;
   
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::Clipped],
    //   p.to_str().unwrap() , false, min_read, min_fail));
    info!( "starting Exon_other" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("ExonOther.tsv");

        run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::ExonOther],  p.to_str().unwrap(), 
                  true, min_read, min_fail, thread as u32, no_beta)?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::ExonOther],
    //   p.to_str().unwrap(), true, min_read, min_fail));
           info!( "starting E_isoform" );
    let mut p = Path::new(&outfile_prefix).to_path_buf();
    let _ = p.set_extension("Isoform.tsv");
            run_one_test( &shared, vec![SplicingCategory::Spliced],
                  vec![SplicingCategory::EIsoform],  p.to_str().unwrap(), 
                  false, min_read, min_fail, thread as u32, no_beta)?;
    //jobs.push((vec![SplicingCategory::Spliced], vec![SplicingCategory::EIsoform],
    //   p.to_str().unwrap() , false, min_read, min_fail));


    //println!("{:?}", );
    let now = Instant::now();
    //let ok = jobs.into_par_iter().
    //try_for_each(|(a, b, c, d, m, mf)| {
    //    let shared_ref = Arc::clone(&shared);
    //    run_one_test(&shared_ref, a, b, c, d, m, mf)
    //});
    //match ok {
    //    Ok(()) => println!("All four tests finished successfully."),
    //    Err(e) => eprintln!("A test failed: {}", e),
    //}

    let elapsed_time = now.elapsed();
    println!("Running slow_function() took {} seconds.", elapsed_time.as_secs());   

        }
    }

    Ok(())
}



       
