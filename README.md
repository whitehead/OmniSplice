Welcome to the omnisplice read me.

# patch note
> 0.5
> rework of compare conditions
> A - Omnisplice statistical approaches have been reworked, and know use R aod library.
you new reauired R with aod instaled in your path.  by default OmniSplice uses  (Fischer extact test, logistic glm). to use beta binomial t-test uses the  "--beta-bin"flag

> B - compare_conditions now as a new parameter '--thread' allowing you to use multi threading. This is usefull because the new test (beta binomial) is way slower. 

> C - I changed the output of "compare_conditions" adding CI, cohens'h, psi-delta, and the position are now represented as chr:start-end(strand). 
> typo Donnor -> Donor corrected. 
> now can selected to 


Omnisplice is a tool to categorize all reads present at exon's extremity.
get more info reading the paper: https://www.biorxiv.org/content/10.1101/2025.04.06.647416v1

This file follow the following organisation:

    - Quick uses
    - 
    - Installation
    - Input Files
    - Command details
    - Output Files
    - Context and discussion.
    - Algorithm overview

Don't hesitate to open an issue if you have any problem or have any suggestions.

## To keep in mind
The same read will be counted as many time as a read overlap different feature.
even it is the same position (different transcript same gene).


## Quick uses
```

the standard pipe-line that I uses:


you first need to run OmniSplice for all your sorted index bam file:
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --anchor 3

then run compares conditions
compare_conditions run-all --control-files control1.junction control2.junction --treatment-files treatment1.junction treatment2.junction  --min-read 10  --outfile-prefix <out file>  --min-read  30 --min-fail 10 --thread 15



# specify the lib type for unstranded data
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --Libtype Unstranded

# extending the overhang to 5 pb
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --overhang 5 


# to get exon_others and associated bam files:
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --read-to-write-junc  exon-other clipped


# for backsplicing you need to first extract the clipped reads
omni_splice --gtf <gtf file> --input <indexBamFile> -o <outFilePrefix> --read-to-write soft-clipped

# then you need bowtie 2 in your OS path, plus a bowtie 2 reference for your genome
backsplicing -i <omnisplice_out.clipped> -o <outputPrefix>> -b <bowtie2 ref> -g <gtf> -m <min clipped size, default 20> 

# for comparison 
compare_conditions run-all --control-files control1.junction control2.junction --treatment-files treatment1.junction treatment2.junction  --min-read 10  --outfile-prefix <out file> ---min-read  30 --min-fail 10 --thread 15

# specific splicing category
compare_conditions run --control-files control1.junction control2.junction --treatment-files treatment1.junction treatment2.junction  --min-read 10  --outfile-prefix <out file> --splicing-fail  Unspliced Clipped Exon_other --splicing-ok Spliced --min-read  30 --min-fail 10 --thread 15


# quick plot of genes (barplot as in manuscript), see "Gene visualisation" below
python gene_visualisation.py --group1 control_S1.junctions control_S2.junctions --group1_name Control --group2 mutant_S1.junctions mutant_S2.junctions --group2_name Mutant --gene_list gene_id1 gene_id2 --outfile_prefix outname
```


## Installation

you need R and Rscript on your path with those library available:
>aod
>data.table
>parallel

to use omni splice you need access to a linux machine (mac never tested).
you need rust, to install rust go to (https://www.rust-lang.org/tools/install). 
you also need cmake. (brew install cmake on mac.)

```
# Clone this repository using git
> git clone <####>
> cd omnisplice
> cargo build --release
# this is a compilation command it will output a bunch of text including some wanring
# and should finish with: "Finished `release` profile [optimized] target(s) "
you can add target directory to your path
```

omnisplice executables will then be at <path>/omnisplice/target/release/omni_splice
if you want to add it to your Path, you can add the release directory to your ~/.bashrc.
example:
export PATH=${PATH}:"<PATH TO>/omnisplice/target/release/"


if you encouter any error, have any questions, want to propose improvment
or let us now by opening a new issue on this github.




### Rust installation
    To see if you have rust installed type "cargo" on the command line, if the terminal return a line with "command not found" you need to install it.
    It is very easy just follow this link: 
        https://www.rust-lang.org/tools/install
    
    if you have rust installed you may need to update it:   
        rustup update


## Input Files
    - gtf: The gtf must have a valid "gene_id" and "transcript_id" for every feature annotated as "exon".
    omnisplice will look at all the exon in this file. if you want to limit your search to subset of gene or exon.
    only include the one your are interested in. (this is particularly helpfull when using the "readToWrite" option as to limit the size of the output).
    - input: a valid position sorted indexed **bam file**


## Command details
    required:
        --input, -i <sortedIndexedBamFile>
        --output-file-prefix, -o <outFilePrefix>
        --gtf, -g <gtfFile>
    optional:
        --libtype <RNA-seq library type, default: frFirstStrand>
        --anchor <minimum anchor on both exon and intron side, default: 3>
        --anchor-exon <minimum anchor on the exon side, overrides --anchor>
        --anchor-intron <minimum anchor on the intron side, overrides --anchor>
        --mapq <minimum MAPQ, default: 13>
        --flag-in <BAM flags a read must have, default: 0>
        --flag-out <BAM flags a read must not have, default: 3840>
        --read-to-write-junc <output reads of specific splicing event categories>
        --read-to-write <(mostly deprecated) output reads by exon-end category>
        --unspliced-def <categories counted as unspliced for splicing efficiency, default: unspliced>
        --spliced-def <categories counted as spliced for splicing efficiency, default: spliced>
        --log-level <default: info>

```
❯ omni_splice --help
Usage: omni_splice [OPTIONS] --input <INPUT> --output-file-prefix <OUTPUT_FILE_PREFIX> --gtf <GTF>

Options:
  -h, --help
          Print help (see a summary with '-h')
  -V, --version
          Print version

Input Options:
  -i, --input <INPUT>
          Name of the input BAM file (sorted and indexed).
  -g, --gtf <GTF>
          Name of the GTF input file defining the features to look at. (v1) Only features
          annotated as exon are considered. If you use --read-to-write or --read-to-write-junc
          with the whole genome, the output can be very large; you may want to subset the
          genes/features you are interested in.
      --libtype <LIBTYPE>
          Library type used for the RNA-seq. Most modern stranded RNA-seq libraries are
          frFirstStrand, which is the default value. Accepted values: frFirstStrand,
          frSecondStrand, fFirstStrand, fSecondStrand, ffFirstStrand, ffSecondStrand,
          rfFirstStrand, rfSecondStrand, rFirstStrand, rSecondStrand, Unstranded,
          PairedUnstranded [default: frFirstStrand]

Output Options:
  -o, --output-file-prefix <OUTPUT_FILE_PREFIX>
          Prefix to be used for the output files.
      --read-to-write <READ_TO_WRITE>...
          Mostly deprecated, use --read-to-write-junc instead. Its only remaining use is
          backsplicing detection with the soft-clipped value. Outputs reads based on their
          exon-end category.
          [possible values: read-through, read-junction, unexpected, fail-pos-filter,
          wrong-strand, fail-qc, empty-pileup, skipped, skipped-unrelated, soft-clipped,
          overhang-fail, empty, all]
      --read-to-write-junc <READ_TO_WRITE_JUNC>...
          Space-separated list of the splicing event categories (relative to the junctions)
          whose reads you want to extract, e.g. spliced, exon-other...
          [possible values: exon-other, spliced, isoform, unspliced, clipped, wrong-strand,
          skipped, skipped-unrelated]
      --unspliced-def <UNSPLICED_DEF>...
          Space-separated list of the categories counted as "unspliced" in the splicing
          efficiency table (this table can be regenerated with the splicing_efficiency
          executable). Accepted values: spliced, unspliced, clipped, exon_other, skipped,
          wrong_strand, isoform. By default, only unspliced (read-through) reads are used.
          To use unspliced and clipped reads: --unspliced-def unspliced clipped
          [default: unspliced]
      --spliced-def <SPLICED_DEF>...
          Space-separated list of the categories counted as "spliced" in the splicing
          efficiency table. Accepted values: spliced, unspliced, clipped, exon_other, skipped,
          wrong_strand, isoform. By default, only spliced reads are used.
          To use spliced and isoform reads: --spliced-def spliced isoform
          [default: spliced]
      --log-level <LOG_LEVEL>
          Log level. Accepted values: error, warn, info, debug, trace [default: info]

Anchor Options:
      --anchor <ANCHOR>
          OmniSplice first identifies reads that align to exon ends and requires them to align
          at least X contiguous bases before the exon boundary. Then, for reads identified as
          unspliced, it requires at least X contiguous bases extending from the exon end into
          the intron. Default for both: 3. Use --anchor to set both values, or set them
          individually with --anchor-exon and --anchor-intron. Must be strictly > 0.
      --anchor-exon <ANCHOR_EXON>
          Minimum number of contiguous aligned bases before the exon boundary.
          Overrides --anchor for the exon side.
      --anchor-intron <ANCHOR_INTRON>
          Minimum number of contiguous aligned bases extending into the intron for unspliced
          reads. Overrides --anchor for the intron side.

QC Options:
      --flag-in <FLAG_IN>
          BAM flags a read must have [default: 0]
      --flag-out <FLAG_OUT>
          BAM flags a read must not have [default: 3840]
      --mapq <MAPQ>
          Minimum mapping quality [default: 13]
```

## Strandness and libtype

Omnisplice by default expect reverse pair end reads.

the libtype argument can be confusing. But it represent your librairy layout and if relevant is helpfull to determine the strandedness of your reads.

- **fr** imply that your paired reads are oriented towards each other (most/majority of modern layout uses this)

    --> <--

- **ff**

    --> -->

- **rf**

    <-- -->


- **FirstStrand** imply that your sequencing reads ar reverse complemented compare to the reference genome

- **SecondStrand** imply that your sequencing reads are NOT reverse complemented compared to the reference genome

### In case of **unstranded** data uses:
- Unstranded (for single end layout)
- PairedUnstranded (for pair end)

OmniSplice uses the flag from the alignment + the layout information (libtype argument) to determine the strandedness of a reads. It does that using a library I made you can check up the code here: https://github.com/rLannes/BAMstrandSpecifier



## Context and discussion.
    OmniSplice does not consume much ressource. Execution time will depend on the size of the bam, and the number of chromosome in the gtf. Memory  will depend on the number of "exon" Feature in the gtf .



## Output Files
    category file (.tsv file): 
        Chromosome | Position of the feature | gene ID | transcript ID | Strand | ExonType | ReadCategory | Read Support
        for each feature they are as many line as ReadCategory found for this specific feature

    with:
        ExonType:   [Donor | Acceptor] indicate if the exon extremity is a junction Donor of acceptor.
            Note: the start of the gene is Acceptor, the end of the gene is Donor.
        
        ReadCategory:
            * ReadJunction(n, m) -> junction read from n, m
            * Unspliced -> read is unspliced
            * Skipped(n, m) -> skipped , read has junction from  n to m that include the feature. 
            * SoftClipped -> Read is clipped at the junction.
            * Empty -> not read found at this junction
            * WrongStrand -> The read come from a fragment which strand is different feom the feature
            * FailPosFilter -> the way the algorithm works to reduce the number of comparison may give some false positive. 
                you can discard it, no meaningfull biologicaly.
            * OverhangFail -> The read fail the overhang test.
            * Unexpected -> if you see this please open a bug report.

    both exon and junction file describe events. But the exon file give the detail by exon end whereas the junction file give the results per junction.

    the exon file (.tsv file):
        contig | gene_name | transcript_name | exon_number | ambiguous | strand | pos | next | exon_type | spliced | unspliced | clipped | exon_other | skipped  | wrong_strand | e_isoform
    
    the junction file (.tsv file):
        contig | gene_name | transcript_name | intron_number | strand | ambiguous | Donor | Acceptor | spliced | unspliced | clipped | exon_other | skipped  | wrong_strand | e_isoform



## Algorithm overview

```
For each read find all exon extremity it does overlap.
for each read R:
    for each extremity E:
        Test if fail pos filter.
        Test if the strand match. -> if False Wrong Strand
        Test if the read skipped the junction.
        Tesf if the read fail the overhang.
        Test if the read is unspliced at the feature specifically.
        Test if the read is SoftClipped.
        Test if the read is a Junction read.
```

## add-on

### Gene visualisation

`gene_visualisation.py` plots, for each transcript of the selected genes, a horizontal stacked barplot
with one bar per intron showing the proportion of reads in each splicing category (as in the manuscript).
It takes the `.junctions` files produced by `omni_splice` and accepts one to three groups (genotypes/conditions);
files from the same group are merged and their counts summed. Groups are plotted side by side.

Requirements: python3 with `numpy`, `pandas` and `matplotlib`. The script imports `common_python` from this repository,
so run it from the repository directory (or add the repository to your `PYTHONPATH`).

```
# one group
python gene_visualisation.py --group1 control_S1.junctions control_S2.junctions --group1_name Control --gene_list gene_id1 --outfile_prefix outname

# two groups, several genes, png output
python gene_visualisation.py --group1 control_*.junctions --group1_name Control --group2 mutant_*.junctions --group2_name Mutant --gene_list gene_id1 gene_id2 --outfile_prefix outname --format .png

# three groups, specific transcripts, custom categories
python gene_visualisation.py --group1 wt_*.junctions --group1_name WT --group2 mutA_*.junctions --group2_name MutA --group3 mutB_*.junctions --group3_name MutB --transcript_list transcript_id1 --splicing_defect SPLICED UNSPLICED CLIPPED EXON_OTHER SKIPPED --outfile_prefix outname
```

One file is written per transcript: `<outfile_prefix>_<gene_id>_<transcript_id><format>`.

    required:
        --group1 <junction files for group 1 (space-separated or glob)>
        --outfile_prefix <prefix for the output files>
        exactly one of:
            --gene_list <gene IDs to plot (space-separated)>
            --transcript_list <transcript IDs to plot (space-separated)>
    optional:
        --group1_name <name for group 1, default: group1>
        --group2 <junction files for group 2>
        --group2_name <name for group 2, default: group2>
        --group3 <junction files for group 3, requires --group1 and --group2>
        --group3_name <name for group 3, default: group3>
        --splicing_defect <categories to plot, default: SPLICED UNSPLICED CLIPPED EXON_OTHER;
                           choices: SPLICED, UNSPLICED, CLIPPED, EXON_OTHER, SKIPPED, SKIPPED_UNR, WRONG_STRAND, E_ISOFORM>
        --color <colors for each category, in the same order as --splicing_defect;
                 must provide at least as many colors as categories>
        --hide_count <do not show the total read count on each bar>
        --format <output format, default: .pdf>
        --fig_width <figure width in inches, default: 6.4>
        --fig_height <figure height in inches, default: 4.8>
        --logging_level, -v <DEBUG, INFO or ERROR, default: INFO>




### Splicing efficiency

`omni_splice` already writes a splicing efficiency table (`<outFilePrefix>.se`) using the `--spliced-def` and `--unspliced-def` options.
The `splicing_efficiency` executable lets you regenerate this table from an existing `.junctions` file with different
definitions of "spliced" and "unspliced", without re-running `omni_splice` on the BAM file.
It is built together with `omni_splice` (`cargo build --release`) and lives in `target/release/`.

```
# default: spliced / (spliced + unspliced)
splicing_efficiency -i <omnisplice_out.junctions> -o <out.se>

# count unspliced and clipped reads as unspliced, spliced and isoform reads as spliced
splicing_efficiency -i <omnisplice_out.junctions> -o <out.se> --unspliced-def unspliced clipped --spliced-def spliced isoform

# also include ambiguous junctions
splicing_efficiency -i <omnisplice_out.junctions> -o <out.se> --ambiguous
```

    required:
        --input, -i <junction file produced by omni_splice>
        --output, -o <output file, must not already exist>
    optional:
        --spliced-def <categories counted as spliced, default: spliced>
        --unspliced-def <categories counted as unspliced, default: unspliced>
        --ambiguous <include ambiguous junctions (exon-intron junctions overlapping another exon); excluded by default>

    accepted categories: spliced, unspliced, clipped, exon_other, skipped, skipped_unrelated, wrong_strand, isoform

The output starts with two comment lines (`#`) recording the spliced and unspliced definitions used, followed by a tab-separated table
sorted by contig, gene and position:

    Contig  start  end  strand  geneID  Ratio  spliced  Unspliced  transcriptID  intronN

- `start` / `end`: donor and acceptor positions of the intron.
- `Ratio`: spliced / (spliced + unspliced), or `NA` if there are no reads.
- `spliced` / `Unspliced`: summed counts of the selected categories.

Junctions without a donor or acceptor position (`.`) are skipped. A temporary `<output>.presorted` file is created and removed at the end.