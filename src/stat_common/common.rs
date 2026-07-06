
use std::hash::{Hash, Hasher};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::{fs, result};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::io::prelude::*;
use std::convert::From;
use std::fmt;
use std::convert::TryFrom;

use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};


use crate::common::utils::ReadAssign;

use super::super::common::error::OmniError;
use super::errors::LogisticRegressionError;
use statrs::distribution::{Discrete, DiscreteCDF, Hypergeometric, StudentsT, ContinuousCDF};
use statrs::statistics::{Min, Max};
use adjustp::{adjust, Procedure};


pub fn apply_bh_correction<T>(
    items: &mut [T],
    get_p: impl Fn(&T) -> Option<f64>,
    set_q: impl Fn(&mut T, f64),
) {
    let mut refs: Vec<&mut T> = Vec::new();
    let mut pvals: Vec<f64> = Vec::new();

    for item in items.iter_mut() {
        if let Some(p) = get_p(item) {
            pvals.push(p);
            refs.push(item);
        }
    }

    if pvals.is_empty() {
        return;
    }

    let qvals = adjust(&pvals, Procedure::BenjaminiHochberg);

    for (r, q) in refs.into_iter().zip(qvals) {
        set_q(r, q);
    }
}



#[derive(Debug, Clone, Copy)]
pub struct TtestResult{
    pub t_stat: Option<f64>,
    pub p_value: Option<f64>,
    pub q_value: Option<f64>
}

impl TtestResult{
    pub fn new_empty() -> Self{
        Self{
            t_stat: None,
            p_value: None,
            q_value: None,

        }
    }
}
pub fn welch_t_test(a: &[f64], b: &[f64]) -> Result<TtestResult, LogisticRegressionError> {

    if a.len() < 2 || b.len() < 2 {
        return Err(LogisticRegressionError::TtestError);
    }
    let mean_a = a.iter().sum::<f64>() / a.len() as f64;
    let mean_b = b.iter().sum::<f64>() / b.len() as f64;
    let var_a = a.iter().map(|x| (x - mean_a).powi(2)).sum::<f64>() / (a.len() as f64 - 1.0);
    let var_b = b.iter().map(|x| (x - mean_b).powi(2)).sum::<f64>() / (b.len() as f64 - 1.0);

    let se = (var_a / a.len() as f64 + var_b / b.len() as f64).sqrt();
    let t_stat = (mean_a - mean_b) / se;

    // Welch-Satterthwaite degrees of freedom
    let df = (var_a / a.len() as f64 + var_b / b.len() as f64).powi(2)
        / ((var_a / a.len() as f64).powi(2) / (a.len() as f64 - 1.0)
            + (var_b / b.len() as f64).powi(2) / (b.len() as f64 - 1.0));

    let t_dist = StudentsT::new(0.0, 1.0, df)?;
    let p_value = 2.0 * (1.0 - t_dist.cdf(t_stat.abs()));

    if !t_stat.is_finite() || !p_value.is_finite() {
        return Err(LogisticRegressionError::TtestError);
    }

    Ok(TtestResult{t_stat: Some(t_stat), p_value: Some(p_value), q_value: None})
}


#[cfg(test)]
mod tttest {
    use super::*;

    #[test]
    fn my_two_group_test() {
        let control_props = vec![0.010, 0.012, 0.011, 0.010, 0.017, 0.010]; // per-replicate proportions
        let treatment_props = vec![0.57, 0.61, 0.56, 0.61, 0.55, 0.58];
        let e = welch_t_test(&control_props, &treatment_props);
        println!("{:?}", e);
    }
}




/// use hyper geom test to compute 2 tailed p-value 
/// the fischer test was returning some negative pvalue?
/// hopefully this works!

pub fn hyper_geom_test(a_succ: u64, a_fail: u64, b_succ: u64, b_fail: u64) -> Result<f64, OmniError>{

    let cond_a = a_succ + a_fail;
    let cond_b = b_succ + b_fail;
    // missing data no meaning
    if (cond_a == 0) && (cond_b == 0){
        return Ok(-1.);
    }
    let succes = a_succ + b_succ;
    let fail = b_fail + a_fail;
    // all succes / all failure no point
    if (fail == 0) && (succes == 0){
        return Ok(1.)
    }
    
    let pop = succes + fail;

    let n = Hypergeometric::new(pop, succes, cond_a).unwrap();
    let k = a_succ;

    let p_k = n.pmf(k);
    let mut p_value: f64 = 0.;
    let log_p_k = n.ln_pmf(k);

    let log_term:Vec<f64> = (n.min()..=n.max())
            .map(|i| n.ln_pmf(i))
            .filter(|&lp| lp <= log_p_k + 1e-10)
            .collect();

    
    if log_term.is_empty(){
        return Ok(1.0)
    }
    let max_lp = log_term.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let p_val = (max_lp + log_term.iter().map(|&lp| (lp - max_lp).exp()).sum::<f64>().ln()).exp();
    return  Ok(p_val.min(1.0))
}




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


pub trait Tester{

    fn to_contengency(&self) -> Vec<u64>{
        let mut res = vec![0 as u64; 4];
        for (i, e) in self.groups().iter().enumerate(){
            match e{
                Genotype::CONTROL => {
                    res[0] += self.success()[i]  as u64 ;
                    res[1] += self.failures()[i] as u64 ;
                },
                Genotype::TREATMENT => {    
                    res[2] += self.success()[i]  as u64 ;
                    res[3] +=  self.failures()[i]  as u64 ;

                }
            }
        }
        res
    }


    fn extract_counts(&self, samples: &Vec<CountsStats>, cat: &Vec<SplicingCategory>, result: &mut Vec<u32>) -> (){
        for count in samples{
            result.push(count.extract_(cat));
        }
    }


    fn format_data(&mut self, treatment: &Vec<CountsStats>, control: &Vec<CountsStats>,
                     successes_cat: &Vec<SplicingCategory>,
                      failures_cat: &Vec<SplicingCategory>) -> (){
        
        let mut success: Vec<u32> = Vec::new();
        let mut failures: Vec<u32>  = Vec::new();
        let mut groups: Vec<Genotype> = Vec::new();
        for _ in 0..control.len(){
            groups.push(Genotype::CONTROL);
        }
        for _ in 0..treatment.len(){
            groups.push(Genotype::TREATMENT);
        }

        self.extract_counts(control, successes_cat, &mut success);
        self.extract_counts(treatment, successes_cat, &mut success);
        self.extract_counts(control, failures_cat, &mut failures);
        self.extract_counts(treatment, failures_cat, &mut failures);

        *self.success_mut() = success;
        *self.failures_mut() = failures;
        *self.groups_mut() = groups;
    }

    fn get_proportion(&self) -> (u32, u32, u32, u32){

        
        let mut ctrl_suc = 0;
        let mut ctrl_fail = 0;
        let mut treat_suc = 0;
        let mut treat_fail = 0;

        for i in 0..self.groups().len(){
            match  self.groups()[i]{
                Genotype::TREATMENT => {          
                        treat_suc += self.success()[i];
                        treat_fail += self.failures()[i];
                    },
                Genotype::CONTROL => {
                        ctrl_suc += self.success()[i];
                        ctrl_fail += self.failures()[i];
                }
            } 

        }
        (ctrl_suc, ctrl_fail, treat_suc, treat_fail)
    }

    fn get_proportion_string(&self) -> (String, String, String, String){

        
        let mut ctrl_suc: Vec<String> = Vec::new();
        let mut ctrl_fail: Vec<String> = Vec::new();
        let mut treat_suc: Vec<String> = Vec::new();
        let mut treat_fail: Vec<String> = Vec::new();

        for i in 0..self.groups().len(){
            match self.groups()[i] {
                Genotype::TREATMENT => {
                    treat_suc.push(self.success()[i].to_string());
                    treat_fail.push(self.failures()[i].to_string());
                }
                Genotype::CONTROL => {
                    ctrl_suc.push(self.success()[i].to_string());
                    ctrl_fail.push(self.failures()[i].to_string());
                }
            }
        }
        (ctrl_suc.join(","), ctrl_fail.join(","), treat_suc.join(","), treat_fail.join(","))
    }

    fn test(&self, donotrun: bool, min_coverage: u32, min_failure: u32) -> TestResults; // donotrun in case we just need to recover prop data but not run test i.e. ambiguous sample
    fn success(&self) -> &Vec<u32>;
    fn success_mut(&mut self) -> &mut Vec<u32>;
    fn failures(&self) -> &Vec<u32>;
    fn failures_mut(&mut self) -> &mut Vec<u32>;
    fn groups(&self) -> &Vec<Genotype>;
    fn groups_mut(&mut self) -> &mut Vec<Genotype>;
}



#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Clone)]
pub enum TestStatus {
    Ok,
    QuasiPerfectSeparation,
    DimensionMistmatch,
    InvalidData,
    InsufficientObservation,
    ambiguous,
    EmptyData,
    ControlIsNull,
    TreatmentIsNull,
    ConvergenceFailed,
    NumericalInstability,
    SingularMatrix,
    PerfectSeparation,
    FisherFallBack,
    TtestFallback,
    HyperGeom,
    FailFilter,
    CIUnavail,
    OddRatioUnavail
}


impl fmt::Display for TestStatus {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {

        match self {
            
            TestStatus::Ok => { write!(f, "Ok") },
            TestStatus::QuasiPerfectSeparation => { write!(f, "QuasiPerfectSeparation") },
            TestStatus::DimensionMistmatch => { write!(f, "DimensionMistmatch") },
            TestStatus::InvalidData => { write!(f, "InvalidData") },
            TestStatus::InsufficientObservation => { write!(f, "InsufficientObservation") },
            TestStatus::ambiguous => { write!(f, "ambiguous") },
            TestStatus::EmptyData => { write!(f, "EmptyData") },
            TestStatus::ControlIsNull => { write!(f, "ControlIsNull") },
            

            TestStatus::HyperGeom =>  { write!(f, "HyperGeomError") },
            TestStatus::TreatmentIsNull => { write!(f, "TreatmentIsNull") },

            TestStatus::ConvergenceFailed => { write!(f, "ConvergenceFailed") },
            TestStatus::NumericalInstability => { write!(f, "NumericalInstability") },
            TestStatus::SingularMatrix => { write!(f, "SingularMatrix") },
            TestStatus::PerfectSeparation => { write!(f, "PerfectSeparation") },
            TestStatus::FisherFallBack => { write!(f, "FisherTest")},
            TestStatus::FailFilter => { write!(f, "Failfilter") }
            TestStatus::CIUnavail => { write!(f, "CIUnavailable") },
            TestStatus::OddRatioUnavail => { write!(f, "ORUnavailable") },
            TestStatus::TtestFallback => { write!(f, "TtestFallback") },

        }
    }
}


impl From<LogisticRegressionError> for TestStatus {
    fn from(item: LogisticRegressionError) -> Self {
        match item{

            LogisticRegressionError::HyperGeomError => TestStatus::HyperGeom,

            LogisticRegressionError::PerfectSeparation{message: _} => TestStatus::PerfectSeparation,
            
            LogisticRegressionError::InvalidProbability(f64) => TestStatus::NumericalInstability ,
            
            LogisticRegressionError::ConvergenceFailure  { iterations: _, final_norm: _ } => TestStatus::ConvergenceFailed,
        
            LogisticRegressionError::SingularMatrix(String) => TestStatus::SingularMatrix,

            LogisticRegressionError::DimensionMismatch { expected: _, got: _ } => TestStatus::DimensionMistmatch,
            
            LogisticRegressionError::InvalidData(String) => TestStatus::InvalidData,
            
            LogisticRegressionError::NumericalInstability(String) => TestStatus::NumericalInstability,

            LogisticRegressionError::EmptyData => TestStatus::EmptyData,

            LogisticRegressionError::ControlIsNull => TestStatus::ControlIsNull,

            LogisticRegressionError::TreatmentIsNull => TestStatus::TreatmentIsNull,

            LogisticRegressionError::CIUnavail(_) => TestStatus::CIUnavail,
            LogisticRegressionError::oddRatioError => TestStatus::OddRatioUnavail,
            LogisticRegressionError::FailUseFisher => TestStatus::FisherFallBack,
            LogisticRegressionError::TtestError => TestStatus::NumericalInstability,
            LogisticRegressionError::TtestStudents(e) => TestStatus::NumericalInstability,
        }
        
    }
}

#[derive(Debug, Clone)]
pub struct TestResults{

    pub control_success: u32,
    pub control_failure: u32,
    pub control_prop: Option<f32>,
    pub predicted_control_prop: Option<f32>,
    pub treatment_success: u32,
    pub treatment_failure: u32,
    pub treatment_prop: Option<f32>,
    pub predicted_treatment_prop: Option<f32>,
    pub status: Option<TestStatus>,
    // Only available if model fit succeeded
    pub p_value: Option<f64>,
    pub q_value: Option<f64>,
    pub odd_ratio: Option<f64>,
    pub or_ci_lower: Option<f64>,
    pub or_ci_upper: Option<f64>,

    pub string_count: (String, String, String, String)
}

impl TestResults{
    pub fn get_empty() -> Self{
        TestResults { control_success: 0, control_failure: 0, control_prop: None,
                      treatment_success: 0, treatment_failure: 0, treatment_prop: None,
                      status: None, p_value: None, q_value: None,
                      predicted_treatment_prop: None, predicted_control_prop: None,
                      odd_ratio: None, or_ci_lower: None,
                      or_ci_upper: None, string_count: ("".to_string(), "".to_string(), "".to_string(), "".to_string())}
    }

    fn helper_(value: Option<f32>) -> String{
        match value{
            Some(v) => format!("{:e}", v),
            None => "nan".to_string()
        }
    }
    fn helper_6(value: Option<f64>) -> String{
        match value{
            Some(v) => format!("{:e}", v),
            None => "nan".to_string()
        }
    }
    fn helper_t(value: &Option<TestStatus>) -> String{
        match value{
            Some(v) => v.to_string(),
            None => "nan".to_string()
        }
    }
    
    pub fn dump_stats(&self, q_val: Option<f64>) -> Vec<String>{

        //format!("{}\t{}\t{}\t{}", 
         vec![           Self::helper_6(self.odd_ratio),
                    Self::helper_6(self.p_value),
                    Self::helper_6(q_val),
                    Self::helper_t(&self.status)]
        //        )
    }
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checking_parse_js() {
        ()
        // I have some weird behavious just checkincg a junction is actulaly picked up by this f()
        /*let mut res : HashMap<String, JunctionStats> = HashMap::new();
        // 
        //
        //
        let control_file = vec!["/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Testv04/OmniSpliceRun_TDP43_anchor3/SRR22002170_R1_001..junctions".to_string(),
        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Testv04/OmniSpliceRun_TDP43_anchor3/SRR22002171_R1_001..junctions".to_string(),
        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Testv04/OmniSpliceRun_TDP43_anchor3/SRR22002172_R1_001..junctions".to_string()];
        let treatment_file = vec!["/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Testv04/OmniSpliceRun_TDP43_anchor3/SRR22002156_R1_001..junctions".to_string(),
        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Testv04/OmniSpliceRun_TDP43_anchor3/SRR22002157_R1_001..junctions".to_string(),
        "/lab/solexa_yamashita/people/Romain/Projets/OmniSplice/Testv04/OmniSpliceRun_TDP43_anchor3/SRR22002158_R1_001..junctions".to_string()];

        for file in control_file{
            info!("parsing {:?}", file);
            parse_js_file(&file, &mut res, Genotype::CONTROL).unwrap();
            info!("done reading");
        }

        for file in treatment_file{
                info!("parsing {:?}", file);
                parse_js_file(&file, &mut res, Genotype::TREATMENT).unwrap();
                info!("done reading");
        }

        let key = format!("{} {} {} {}", "7", "+", "127173138",	"127176706");
        println!("{:?}", res.get(&key));

        let j = res.get(&key).unwrap();
        //for (k, j) in junction.iter(){
        let mut glm = GLM::new(  &j.control_count,
                                        &j.treat_count,
                                                &successes_cat,
                                                &failures_cat,
                                            k.to_owned());
        x = glm.test(true);
        
    }*/}

//key = format!("{} {} {} {}", contig, strand, start, end);
}