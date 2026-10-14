#![allow(unused)]
use clap::CommandFactory;
use clap::Parser;

use CigarParser::cigar::Cigar;
use rust_htslib::bam::header;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::{IndexedReader, Read};
use std::collections::HashMap;
use std::collections::HashSet;
use std::default;
use std::error::Error;
use std::fmt::format;
use std::fs::File;
use std::hash::Hash;
use std::io::BufWriter;
use std::io::prelude::*;
use strand_specifier_lib::{LibType, check_flag};

use std::ffi::OsStr;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
mod common;
use common::error::OmniError;
use common::it_intron::dump_tree_into_raw_exon_junction;
use std::convert::From;
use std::fs;
use std::str;

use bio::data_structures::interval_tree::IntervalTree;

use flexi_logger::{FileSpec, Logger, WriteMode};
use log::{debug, error, info, trace, warn};

//use crate::common::utils::ReadAssign;
use crate::common::gtf_::{
    get_all_junction_for_a_gene, get_invalid_pos, get_junction_from_gtf, gtf_to_hashmap, exon_intervalltree
};
use crate::common::it_intron::TreeDataIntron;
use crate::common::it_intron::{
    dump_tree_to_cat_results, interval_tree_from_gtfmap, update_tree_from_bam,
};
//use crate::common::it_approches::{
//    dump_tree_to_cat_results, gtf_to_tree, update_tree_with_bamfile,
//};
//use crate::common::point::{read_gtf, InsideCounter, PointContainer};
use crate::common::junction_file::junction_file_from_table;
use crate::common::read_record::file_to_table;
use crate::common::utils;
use crate::common::utils::Exon;
use crate::common::utils::SplicingEvent;
use crate::common::utils::{ReadsToWriteSEvent,ReadToWriteHandleJunc,
    ExonType, ReadAssign, ReadToWriteHandle, ReadsToWrite, update_read_to_write_handle, update_read_to_write_handle_junc
};
mod splicing_efficiency;

// TODO add LibType

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Name of the input BAM file (sorted and indexed).
    #[arg(short, long, required = true, help_heading = "Input Options")]
    input: String,
    /// Prefix to be used for the output files.
    #[arg(short, long, required = true, help_heading = "Output Options")]
    output_file_prefix: String,
    /// Name of the GTF input file defining the features to look at.
    /// (v1) Only features annotated as exon are considered.
    /// If you use --read-to-write or --read-to-write-junc with the whole genome, the output can be very large;
    /// you may want to subset the genes/features you are interested in.
    #[arg(short, long, required = true, help_heading = "Input Options")]
    gtf: String,


    /// OmniSplice first identifies reads that align to exon ends and requires them to align at least X contiguous bases before the exon boundary. Then, for reads identified as unspliced, it requires at least X contiguous bases extending from the exon end into the intron.
    /// Default for both: 3.
    /// Use --anchor to set both values, or set them individually with --anchor-exon and --anchor-intron.
    /// Must be strictly > 0.
    #[arg(long, help_heading = "Anchor Options")]
    anchor: Option<i64>,
    /// Minimum number of contiguous aligned bases before the exon boundary.
    /// Overrides --anchor for the exon side.
    #[arg(long, help_heading = "Anchor Options")]
    anchor_exon: Option<i64>,
    /// Minimum number of contiguous aligned bases extending into the intron for unspliced reads.
    /// Overrides --anchor for the intron side.
    #[arg(long, help_heading = "Anchor Options")]
    anchor_intron: Option<i64>,

    /// BAM flags a read must have.
    #[arg(long, default_value_t = 0, help_heading = "QC Options")]
    flag_in: u16,
    /// BAM flags a read must not have.
    #[arg(long, default_value_t = 3840, help_heading = "QC Options")]
    flag_out: u16,
    /// Minimum mapping quality.
    #[arg(long, default_value_t = 13, help_heading = "QC Options")]
    mapq: u8,
    /// Mostly deprecated, use --read-to-write-junc instead. Its only remaining use is backsplicing detection with the soft-clipped value.
    /// Outputs reads based on their exon-end category.
    #[clap(long, value_parser, value_delimiter = ' ', num_args = 1.., help_heading = "Output Options")]
    read_to_write: Vec<ReadsToWrite>,
    /// Space-separated list of the splicing event categories (relative to the junctions) whose reads you want to extract,
    /// e.g. spliced, exon-other...
    #[clap(long, value_parser, value_delimiter = ' ', num_args = 1.., help_heading = "Output Options")]
    read_to_write_junc: Vec<ReadsToWriteSEvent>,
    /// Space-separated list of the categories counted as "unspliced" in the splicing efficiency table
    /// (this table can be regenerated with the splicing_efficiency executable).
    /// Accepted values: spliced, unspliced, clipped, exon_other, skipped, wrong_strand, isoform.
    /// By default, only unspliced (read-through) reads are used.
    /// To use unspliced and clipped reads: --unspliced-def unspliced clipped
    #[clap(long, value_parser, default_value = "unspliced", value_delimiter = ' ', num_args = 1.., help_heading = "Output Options")]
    unspliced_def: Vec<String>,

    /// Space-separated list of the categories counted as "spliced" in the splicing efficiency table.
    /// Accepted values: spliced, unspliced, clipped, exon_other, skipped, wrong_strand, isoform.
    /// By default, only spliced reads are used.
    /// To use spliced and isoform reads: --spliced-def spliced isoform
    #[clap(long, value_parser, default_value = "spliced", value_delimiter = ' ', num_args = 1.., help_heading = "Output Options")]
    spliced_def: Vec<String>,

    /// Library type used for the RNA-seq. Most modern stranded RNA-seq libraries are frFirstStrand, which is the default value.
    /// Accepted values: frFirstStrand, frSecondStrand, fFirstStrand, fSecondStrand, ffFirstStrand, ffSecondStrand, rfFirstStrand,
    /// rfSecondStrand, rFirstStrand, rSecondStrand, Unstranded, PairedUnstranded
    #[clap(long, value_parser, default_value = "frFirstStrand", help_heading = "Input Options")]
    libtype: LibType,

    /// Log level. Accepted values: error, warn, info, debug, trace
    #[clap(long, value_parser, default_value = "info", help_heading = "Output Options")]
    log_level: String,
}

/// Runs the core of the program: parses a GTF and a BAM file, and writes a category file and, if requested, a read file.
fn main_loop(
    //output: String,
    //out_j: String,
    gtf: String,
    bam_input: String,
    anchor_exon: i64,
    anchor_intron: i64,
    flag_in: u16,
    flag_out: u16,
    mapq: u8,
    output_write_read_handle: &mut ReadToWriteHandle,
    output_write_read_handle_jun: &mut ReadToWriteHandleJunc,
    librairy_type: LibType,
    gtf_hashmap: &HashMap<String, HashMap<String, Vec<Exon>>>,
    valid_j_gene: &HashMap<String, HashSet<(i64, i64)>>,
) -> Result<HashMap<String, IntervalTree<i64, TreeDataIntron>>, OmniError> {
    info!("Launching the main loop parsing the bam.");
    let bam_file = bam_input;

    let gtf_file = gtf;
    info!("Building the interval tree.");
    let mut hash_tree =
        interval_tree_from_gtfmap(gtf_hashmap).expect("failed to generate the hash tree from gtf");

    /// TODO duplication here! done that do in MAIN!
    //let junction_ambi = get_invalid_pos(&gtf_file)?;
    let junction_ = get_junction_from_gtf(&gtf_file, &librairy_type)?;


    update_tree_from_bam(
        &mut hash_tree,
        &bam_file,
        librairy_type, //LibType::frFirstStrand,
        anchor_exon,
        anchor_intron,
        flag_in,
        flag_out,
        mapq,
        output_write_read_handle,
        output_write_read_handle_jun,
        &junction_,
        valid_j_gene,
    );

    return Ok(hash_tree);

}

fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
    let args = Args::parse();

    match args.libtype {
        LibType::Invalid => panic!("invalid librairy type"),
        _ => (),
    }
    let mut output_file_prefix = args.output_file_prefix.clone();

    let out_exons = format!("{}{}", output_file_prefix, ".exons");
    let out_raw = format!("{}{}", output_file_prefix, ".raw");
    let splicing_defect = format!("{}{}", output_file_prefix, ".se");
    let out_junction = format!("{}{}", output_file_prefix, ".junctions");
    let log_file = format!("{}{}", output_file_prefix, ".log");

    // anchor option
    let anchor = match args.anchor{
        Some(x) => x,
        None => 3 as i64
    };
    let anchor_exon = match args.anchor_exon{
        Some(x) => x,
        None => anchor,
    };
    let anchor_intron = match  args.anchor_intron{
        Some(x) => x,
        None => anchor,
    };


    let mut clipped = false;

    let intermediate = output_file_prefix.to_owned();
    let m = Path::new(&intermediate);

    let mut log_dir = m.parent().unwrap().to_str().unwrap().to_string();
    if log_dir.is_empty() {
        log_dir = "./".to_string();
    }
    println!("logging to directory : {} ", log_dir);
    Arc::new(
        Logger::try_with_str(&args.log_level)? // set the default log level
            .log_to_file(
                FileSpec::default()
                    .directory(log_dir) //m.parent().unwrap().to_str().unwrap().to_string())          // create files in folder ./log_files
                    .basename(m.file_stem().unwrap().to_str().unwrap().to_string())
                    .discriminant("OmniSplice") // use infix in log file name
                    .suffix("log"),
            )
            .write_mode(WriteMode::Async)
            //.format(my_formatter)                  // optional custom format
            .write_mode(WriteMode::SupportCapture)
            .start()
            .expect("Failed to initialise logger"),
    );
    info!("ARGS:");
    info!("\tInput file: {}", args.input);
    info!("\tOutput prefix: {}", args.output_file_prefix.clone());
    info!("\tGTF file: {}", args.gtf);
    info!("\tAnchor: {}", anchor);
    info!("\tAnchor_exon: {}", anchor_exon);
    info!("\tAnchor_intron: {}", anchor_intron);
    info!("\tFlag in: {}", args.flag_in);
    info!("\tFlag out: {}", args.flag_out);
    info!("\tMAPQ: {}", args.mapq);
    info!("\tReads to write: {:?}", args.read_to_write);
    info!("\tReads to write junction: {:?}", args.read_to_write_junc);
    info!("\tUnspliced def: {:?}", args.unspliced_def);
    info!("\tSpliced def: {:?}", args.spliced_def);
    info!("\tLibrary type: {:?}", args.libtype);
    info!("\tLog level: {}", args.log_level);
    info!("END ARGS:");

    if anchor <= 0 || anchor_exon <= 0 || anchor_intron <= 0 {
        error!("anchor options must be strictly positive (> 0)");
        return Err(OmniError::Expect("Error: anchor options expect strictly positive values (>= 1)".to_string()).into());
    }
    
    let header_reads_handle_exon = "contig\taln_start\tgene_id\ttranscript_id\tstrand\texon_JType\tread_name\tcig\tsequence\n".as_bytes();
    let mut read_out_handle = ReadToWriteHandle::new();
    update_read_to_write_handle(
        &mut read_out_handle,
        args.read_to_write,
        header_reads_handle_exon,
        &output_file_prefix,
    );

    let bam = rust_htslib::bam::Reader::from_path(&args.input).unwrap();
    let header = rust_htslib::bam::Header::from_template(bam.header());
    drop(bam);

    let header_handle_junc = "contig\tgene_id\ttranscript_id\tstrand\tJ_left\tJ_right\tread_name\taln_start\tcigar\tflag\tsequence\n".as_bytes();
    let mut read_out_handle_jun = ReadToWriteHandleJunc::new();
        update_read_to_write_handle_junc(
        &mut read_out_handle_jun,
        args.read_to_write_junc,
        header_handle_junc,
        &output_file_prefix,
        &header
    );


    info!("Parsing gtf file.");
    info!("Getting all ambiguous positions.");
    let ambiguous_position = get_invalid_pos(&args.gtf.clone())?;
    // for ambigious
    let exon_it = exon_intervalltree(&args.gtf.clone())?;
    let gtf_hashmap = gtf_to_hashmap(&args.gtf.clone()).expect("failed to parse gtf");
    info!("Getting all valid junctions to identify isoforms.");
    let valid_j_gene = get_all_junction_for_a_gene(&gtf_hashmap)
        .map_err(|e| format!("failed to get all junctions for genes, must abort: {e}"))?;

    let tree = main_loop(
        // output.clone(),
        //  junction_file.clone(),
        args.gtf.clone(),
        args.input,
        anchor_exon,
        anchor_intron,
        args.flag_in,
        args.flag_out,
        args.mapq,
        &mut read_out_handle,
        &mut read_out_handle_jun,
        args.libtype,
        &gtf_hashmap,
        &valid_j_gene,
    )?;

    info!("Main loop ended, writing results.");

    let junction_order: Vec<SplicingEvent> = vec![
        SplicingEvent::Spliced,
        SplicingEvent::Unspliced,
        SplicingEvent::Clipped,
        SplicingEvent::ExonOther,
        SplicingEvent::Skipped,
        SplicingEvent::SkippedUnrelated,
        SplicingEvent::WrongStrand,
        SplicingEvent::Isoform,
    ];

    dump_tree_into_raw_exon_junction(
        &tree,
        &out_raw,
        &out_exons,
        &out_junction,
        &ambiguous_position,
        &exon_it,
        &junction_order,
    )?;


    //junction_file_from_table(&table, &junction_file);
    splicing_efficiency::to_se_from_junction(
        &out_junction.clone(),
        &splicing_defect,
        args.spliced_def,
        args.unspliced_def,
        false
    );
    Ok(())
}

//// TODO OPTIMIZATION:
//// CHANGE hasmap algorithm
//// refactor gtf
//// use multiThreading (1 tree chromosomes) in //?
