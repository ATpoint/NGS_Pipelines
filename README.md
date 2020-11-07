# NGS_Pipelines

Bash-based pipelines for processing of ATAC/ChIP-seq, RNA-seq and scRNA-seq data.
Not intended for public use, therefore use at your own risk, without any warranty.

## Software Installation
It is recommended to install software via the `miniconda` package manager into a separate environment.
A Linux system and `miniconda3` is assumed, on Mac it should work as well, but not extensively tested.


```
conda create --name Pipelines
conda activate Pipelines

Tools=(bedtools=2.29.2 bioconductor-edger=3.32.0 bowtie2=2.4.2 coreutils cutadapt \
       fastqc subread macs2=2.2.7.1 mawk multiqc r-base=4.0 parallel picard pigz \
       salmon=1.3.0 samblaster samtools=1.11 seqtk tabix ucsc-bedgraphtobigwig)
        
printf '%s\n' "${Tools[@]}" > install_software.txt

conda install --file install_software.txt
```

In order to activate the envir from inside a SLURM script add these lines to the top of the script as 
suggested at [this Github issue](https://github.com/conda/conda/issues/7980).

```
eval "$(conda shell.bash hook)"
conda activate <env-name>
```

## Available Pipeline

<br>
Run any of the bash scripts without arguments or with `-h/--help` to see the help section with all available arguments.
<br>

#### `DNAseq_v1.0.0.sh`

Pipeline for alignment, filtering and QC/FRiP assessment of DNA-seq assays such as ChIP-seq and ATAC-seq.

```{bash}
------------------------------------------------------------------------------------------------------------------
-h | --help          : show this message                                                           {}
-g | --genome        : path to genome index prefix                                                 {}
-j | --jobnumber     : number of parallel jobs (GNU parallel)                                      {4}
-t | --threads       : number of threads per job for alignment                                     {16}          
-f | --format        : type of input files (fq_se/pe, bam_se/pe)                                   {}           
-a | --atacseq       : turn on ATAC-seq mode                                                       {FALSE}
-o | --onlyAln       : perform only the alignment and then exit                                    {FALSE}
-m | --chrM          : name of the mitochondrial chromosome (for ATAC-seq % reads mapping to it)   {chrM}
-l | --memSort       : memory per thread for BAM sorting by samtools in the form of e.g. "1G"      {1G}
-y | --threadsSort   : threads for sorting operations                                              {8}
-c | --chrRegex      : use this regex to select chromosomes to keep during filtering               {chr[1,9,X,Y]}
-k | --keepDup       : do not remove PCR duplicates from alignment                                 {FALSE}  
-q | --minMAPQ       : keep only alignments with that minimal MAPQ                                 {20}
-x | --checktools    : check whether required software is in PATH                                  {}
-d | --noalignment   : do not run alignment                                                        {FALSE}
-s | --nofrips       : do not run the QC assessment via peak calling and FRiP calculation          {FALSE}
-b | --layout        : whether se or pe, only relevant if --noalignment is set (se,pe)             {se}
-w | --fripqcJobs    : number of parallel jobs for the QC/FRiP part                                {36}
------------------------------------------------------------------------------------------------------------------
```

Input data are gzipped fastq or unaligned BAMs with the naming convention:
- single-end: `Basename.fastq.gz`
- paired-end: `Basename_1.fastq.gz`, `Basename_2.fastq.gz`
- ubam se/pe: `Basename_ubam.bam`

Workflow includes adapter trimming, alignment, removal of non-primary, duplicated and low-MAPQ reads as well as alignments on non-primary chromosomes.
Output are the unfiltered sorted alignment (`Basename_raw.bam`) and the filtered one `Basename_dedup.bam`. The user can specify which chromosomes shall be retained during the filtering via `--chrRegex` which is basically the regex used by `grep` on the chromosome names. Default is `chr[1,9,X,Y]` which will keep (here intended for the mouse and human genome) all chromosomes prefixed with chr followed by a number, plus the sec chromosomes X and Y, but discard the unplaced contigs etc.
In ATAC-seq mode (`--atacseq`) it also outputs a BED file with the transposase cutting sites (shifted +4/-5) as both compressed BED and bigwig,
with the names `Basename_cutsites.bed.gz` and `Basename_cutsites.bigwig`.
The pipeline can also perform some basic QC by calling peaks with `macs2` and then calculate the Fraction Of Reads in Peaks (FRiPs) as an estimate of the signal/noise ratio. The FRiPs per sample are then in `FRiPs.txt`. For ATAC-seq there will also be `mtDNA_percent.txt` which contains the percentage of reads per sample mapping to the mitochondrial chromosome (specified via `--chrM`). If one already has alignments from this pipeline one can skip it and only perform the FRiP QC via `--noalignment`, and one can skip the FRiP QC when only running the alignments via `--nofrips`.

As minimum input the path to the `bowtie2` index files must be provided via `--genome` as well as the format (fq_se, fq_pe, bam_se, bam_pe) to indicate input format and sequencing layout (single, paired-end).
Run with `--checktools` to check whether all reqiured software is in `$PATH`. If not `missing_tools.txt` will contain the names of the missing tools.
That check is automatically performed (if not specified explicitely) before every run.

The input files are all fastq.gz or uBAM files in the directory where the script sits. They don't have to be (or can be) explicitely specified.

<br>
<br>

#### `RNAseq_SA_v1.0.0.sh`

The RNA-seq pipeline using `salmon` for quantification of fastq files against a transcriptome.
Run script without arguments to display this help message:

```{bash}
---------------------------------------------------------------------------------------------

-h | --help        : Show this message                               {}
-i | --idx         : the transcriptome index folder                  {}
-m | --mode        : single or paired-end data (single,paired)       {}
-n | --noLength    : turn off length correction for end-tagged
                     libraries                                       {FALSE}
-t | --threads     : number of threads per run                       {16}
-j | --njobs       : number of parallel jobs for salmon              {4}
-l | --libtype     : library type                                    {A}

-s | --fldMean     : mean insert size for single-end data            {250}    
-q | --fldSD       : standard deviation for --fldMean                {25}
-a | --additional  : any additional salmon arguments                 {}    
-c | --trim        : whether to trim adapters via cutadapt           {FALSE}
-d | --adapter     : the adapter sequence to trim, default is TruSeq {AGATCGGAAGAGC}
-y | --trimthreads : threads per job for cutadapt                    {2}
-x | --trimjobs    : GNU parallel jobs for cutadapt                  {10}

---------------------------------------------------------------------------------------------
```

Input files are gzipped fastq, either single-or paired-end with naming concentions as in the DNAseq pipeline above,
so `Basename.fastq.gz` for single-end and `Basename_1.fastq.gz`/`Basename_2.fastq.gz` for paired-end data.
All files with this suffix in the current directory of the script will be used as input.
Script checks whether reqiured tools are in `$PATH` before scarting the job. Missing tools are in `missing_tools.txt`

The pipeline can optionally run `cutadapt` for adapter trimming.
It will always estimate mapping uncertainty using Gibbs sampling (100x).
If using an index that contains the entire human or mouse genome as decoy one should not run more than four parallel 
quantification jobs (`--njobs`) on the standard 92GB HPC nodes. Four jobs with 16 threads each usually works well.
Note that the value goven to `--trimthreads` must be multiplied by two (for single-end) and three (for paired-end) data
as `cutadapt` will pass that parameter to `pigz` for compression of the output files. The defaults would therefore need about 60 cores.

<br>
<br>

#### `scRNAseq_v1.0.0.sh`

The scRNA-seq pipeline for droplet-based data using `alevin` for quantification, CB and UMI extraction/deduplication. 
Run script without arguments to display this help message:

```{bash}
-h | --help       : Show this message                                                                    
-i | --idx        : the transcriptome index                              {}                                                                                                    
-c | --chemistry  : the kit chemistry (dropseq, chromium, chromiumV3)    {chromiumV3}
-g | --tgmap      : transcript2gene conversion table                     {}                                
-m | --mtgenes    : list with mitochondrial genes for CB whitelisting    {}                                    
-r | --rrnagenes  : list with rRNA genes for CB whitelisting             {}                                
-l | --libtype    : library type                                         {ISR}                             
-t | --threads    : number of threads per sample                         {36}
-j | --njobs      : number of parallel jobs                              {2}
-a | --additional : any additional parameters to pass to alevin          {}
```

CBs and UMIs are expected in R1 and cDNA is R2.
No more than two samples should be quantified in parallel when using the standard 92GB HPC nodes.
We usually use an index that harbors both spliced- and unspliced transcripts plus the entire genome as decoy,
therefore memory footprint is considerable. By default no mapping uncertainty via bootstrapping are produced.
If this is desired then feed this in via the `--additional` argument.

<br>
<br>

#### `Bam2Bigwig_v1.0.0.sh`

Accepts bam files as input and produces bigwig files:

```{bash}
-h | --help          : show this message                                                           {}
-b | --bams          : space-delimited string of input bam files in double quotes or "*.bam"       {}
-m | --mode          : <single> or <paired>, if paired will use -pc option in bedtools genomecov   {}
-a | --atacseq       : use +4/-5 shifted 5-prime end of reads to calcualte coverage                {FALSE}
-e | --extend        : numeric value to extend reads to fragment size, see details                 {0}
-n | --normalize     : if set then normalizes bedGraphs using TMM from edgeR, see details          {FALSE}
-u | --useexisting   : use existing scaling_factors.txt to grep SFs from                           {FALSE}
-p | --peaks         : peaks in BED format                                                         {}
-j | --njobs         : number of parallel (GNU parallel) jobs to create bedGraphs, see details.    {1}
-t | --threads       : number of threads for featureCounts (if --normalize)                        {1}
-q | --sortthreads   : number of threads for sorting the bedGraph                                  {1}
-w | --sortmem       : memory for sorting the bedGraph                                             {1G}
```

The workhorse is `bedtools genomecov` to create the genome-wide pileups. Several options to customize results are available.
In ATAC-seq mode (`--atacseq`) only the 5' ends shifted by +4/-5, that is the transposase cutting sites are counted. 
One can extend reads to fragments using `--extend`, or in ATAC-seq mode use the same option to extend the cutting site by the provided 
value in both directions, e.g. 50 to get a total window of 100bp. The latter is done by `bedtools slop`, read extension is simply the `-fs` option of `genomecov`. 
For paired-end data (`--mode paired`) the TLEN is used to connect both mates getting the actual fragment coverage (option `-pc` of `genomecov`).
If one passes a BED file with intervals as `--peaks` together with `--normalize`, then it will produce a count matrix based on the bam files for these intervals and then uses `edgeR` to calculate per-sample scaling factors using the TMM normalization method. 
The resulting factors are then used to divide the score (that is column 4 of a bedGraph) by.
For productive use one can easily import bigwigs into `R` using `rtacklayer::import` which returns a GRanges object 
with the coverage as "score" column. Bigwigs can be visualized in the IGV viewer as well in a memory-efficient fashion (unlike bedGraph).
If one aims to average bigwigs there are two options:
1) Use `wiggletooms mean` which is cumbersome because it returns wig, which needs to bed converted back to bedGraph and then back to bigwig.
If one converts the wig directly to bigwig it will produce a very large files, much larger than bigwig produced from bedGraph, no idea why.
2) Convert bigwig back to bedGraph or bedGraph.gz and then use `bedtools unionbedg`. This could be done in a stream-like fashion like:
`bedtools unionbedg -i <(bigWigToBedGraph in1.bigwig /dec/stdout) <(bigWigToBedGraph inN.bigwig /dec/stdout) (...)`.



