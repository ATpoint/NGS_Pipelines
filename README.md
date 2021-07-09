# NGS_Pipelines

Pipelines for precessing of ATAC/ChIP-seq, RNA-seq and scRNA-seq data.

Note that these pipelines are for internal use and come without any warranty. Default settings are tailored for use on a 72-core HPC node with > 80GB RAM.
it is always expected that the input files (fastq, UBAM/uCRAM) are in `$(pwd)`.

## Software

- the `environment.yml` contains a conda environment (built on CentOS-7) with all required software

- a Docker container based on that environment is available from the [Docker Hub](https://hub.docker.com/r/atpoint/phd_project) based on this [Github repo](https://github.com/ATpoint/phd_project_docker).

The template submission script `submission_scripts/generic_template.slurm` contains instructions on how to launch the workflows via SLURM + Singularity.

## Available Pipeline

<br>
Run any of the bash scripts without arguments or with `-h/--help` to see the help section with all available arguments.
<br>

#### `DNAseq_v1.0.2.sh`

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

Input data are gzipped fastq or unaligned BAMs/CRAMs with the naming convention:  
- single-end: `Basename.fastq.gz`
- paired-end: `Basename_1.fastq.gz`, `Basename_2.fastq.gz`
- ubam se/pe: `Basename_ubam.bam`

Workflow includes adapter trimming, alignment, removal of non-primary, duplicated and low-MAPQ reads as well as alignments on non-primary chromosomes.
Output are the unfiltered sorted alignment (`basename_raw.bam`) and the filtered one `basename_dedup.bam`. The user can specify which chromosomes shall be retained during the filtering via `--chrRegex` which is basically the regex used by `grep` on the chromosome names. Default is `chr[1,9,X,Y]` which will keep (here intended for the mouse and human genome) all chromosomes prefixed with chr followed by a number, plus the sec chromosomes X and Y, but discard the unplaced contigs etc.
In ATAC-seq mode (`--atacseq`) it also outputs a BED file with the transposase cutting sites (shifted +4/-5) as both compressed BED and bigwig,
with the names `Basename_cutsites.bed.gz` and `Basename_cutsites.bigwig`.  
The pipeline can also perform some basic QC by calling peaks with `macs2` and then calculate the Fraction Of Reads in Peaks (FRiPs) as an estimate of the signal/noise ratio. The FRiPs per sample are then in `FRiPs.txt`. For ATAC-seq there will also be `mtDNA_percent.txt` which contains the percentage of reads per sample mapping to the mitochondrial chromosome (specified via `--chrM`). If one already has alignments from this pipeline one can skip it and only perform the FRiP QC via `--noalignment`, and one can skip the FRiP QC when only running the alignments via `--nofrips`.

As minimum input the path to the `bowtie2` index files must be provided via `--genome` as well as the format (fq_se, fq_pe, bam_se, bam_pe) to indicate input format and sequencing layout (single, paired-end).
Run with `--checktools` to check whether all required software is in `$PATH`. If not `missing_tools.txt` will contain the names of the missing tools.
That check is automatically performed (if not specified explicitely) before every run.

The input files (fastq.gz/uBAM/uCRAM) are expected in the same directory as the script.

For submission via SLURM and running via Singularity with our [Docker image](https://hub.docker.com/r/atpoint/phd_project) one could use:

```bash


singularity_basic_pipelines="singularity exec --bind=path/to/dir-to-mount <image.sif> echo '' && ulimit -u 50000"

echo -e '#!/bin/bash'"\n"'eval ${singularity_basic_pipelines} && eval "$(conda shell.bash hook)" && conda activate Pipelines' > submit.slurm \
&& echo 'ls *.gz | parallel -j <paralleljobs> "fastqc {}"' >> submit.slurm \
&& echo "./DNAseq_v1.0.X.sh --genome ${idx_bowtie2}/$(basename ${ref_genome}) --atacseq --format fq_pe" >> submit.slurm

#/ submit to SLURM e.g:
sbatch --nodes=1 --ntasks-per-node=72 --mem=80G --partition=normal --job-name=jobname --time=24:00:00 submit.slurm

```

After running the pipeline one can use the `cleanup.sh` script in this repo to sort output into folders.
  
<br>
<br>

#### `RNAseq_v1.0.2.sh`

The RNA-seq pipeline using `salmon` for quantification of fastq files against a transcriptome.
Run script without arguments to display this help message:

```{bash}
------------------------------------------------------------------------------------------------------------------

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

------------------------------------------------------------------------------------------------------------------
```

Input files are gzipped fastq, either single-or paired-end with naming concentions as in the DNAseq pipeline above,
so `Basename.fastq.gz` for single-end and `Basename_1.fastq.gz`/`Basename_2.fastq.gz` for paired-end data.
All files with this suffix in the current directory of the script will be used as input.
Script checks whether reqiured tools are in `$PATH` before scarting the job. Missing tools are in `missing_tools.txt`

The pipeline can optionally run `cutadapt` for adapter trimming.
If using an index that contains the entire human or mouse genome as decoy one should probably not run more than four parallel 
quantification jobs (`--njobs`) on the standard HPC nodes but this is not extensively tested. 
Four jobs with 16 threads each usually works well without touching memory limits.
Note that the value goven to `--trimthreads` must be multiplied by two (for single-end) and three (for paired-end) data
as `cutadapt` will pass that parameter to `pigz` for compression of the output files. The defaults would therefore need about 60 cores.

<br>
<br>

#### `scRNAseq_v1.0.0.sh`

Deprecated, see https://github.com/ATpoint/sc_preprocess.  

The scRNA-seq pipeline for droplet-based data using `alevin` for quantification, CB and UMI extraction/deduplication. 
Run script without arguments to display this help message:

```{bash}
------------------------------------------------------------------------------------------------------------------
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
------------------------------------------------------------------------------------------------------------------
```

CBs and UMIs are expected in R1 and cDNA is R2.
Not extensively tested, but probably no more than two or three samples should be quantified in parallel when using the standard HPC nodes to avoid memory overload.
We usually use an index that harbors both spliced- and unspliced transcripts plus the entire genome as decoy,
therefore memory footprint is quiet extensive (compared to a normal txtome-only index).
For creation of such an index see: https://github.com/ATpoint/SingleCell
If one wishes to pass additional arguments to alevin this is possible via the `--additional` argument in double quotes, 
e.g `--additional "--numCellGibbsSamples"`.

<br>
<br>

#### `Bam2Bigwig_v1.0.0.sh`

Accepts bam files as input and produces bigwig files:

```{bash}
------------------------------------------------------------------------------------------------------------------
-h | --help          : show this message                                                           {}
-b | --bams          : space-delimited string of input bam files in double quotes or "*.bam"       {}
-m | --mode          : <single> or <paired>, if paired will use -pc option in bedtools genomecov   {}
-a | --atacseq       : use +4/-5 shifted 5-prime end of reads to calcualte coverage                {FALSE}
-e | --extend        : numeric value to extend reads to fragment size, see details                 {0}
-n | --normalize     : if set then normalizes bedGraphs using TMM from edgeR, see details          {FALSE}
-k | --useexistingcm : use this existing count matrix to calculate SFs from                        {FALSE}
-u | --useexistingsf : use existing scaling_factors.txt to grep SFs from                           {FALSE}
-p | --peaks         : peaks in BED format                                                         {}
-j | --njobs         : number of parallel (GNU parallel) jobs to create bedGraphs, see details.    {1}
-t | --threads       : number of threads for featureCounts (if --normalize)                        {1}
-q | --sortthreads   : number of threads for sorting the bedGraph                                  {1}
-w | --sortmem       : memory for sorting the bedGraph                                             {1G}
------------------------------------------------------------------------------------------------------------------
```

The workhorse is `bedtools genomecov` to create the genome-wide pileups. Several options to customize results are available.
In ATAC-seq mode (`--atacseq`) only the 5' ends shifted by +4/-5, that is the transposase cutting sites are counted. 
One can extend reads to fragments using `--extend`, or in ATAC-seq mode use the same option to extend the cutting site by the provided 
value in both directions, e.g. 50 to get a total window of 100bp. The latter is done by `bedtools slop`, read extension is simply the `-fs` option of `genomecov`. 
For paired-end data (`--mode paired`) the TLEN is used to connect both mates getting the actual fragment coverage (option `-pc` of `genomecov`).
If one passes a BED file with intervals as `--peaks` together with `--normalize`, then it will produce a count matrix based on the bam files for these intervals and then uses `edgeR` to calculate per-sample scaling factors using the TMM normalization method. 
The resulting factors are then used to divide the score (that is column 4 of a bedGraph) by.
Function can also use existing factors (`--useexistingsf`) which then have to be in a file `scaling_factors.txt` or use an existing count matrix
via `--useexistingcm`. 
For productive use one can easily import bigwigs into `R` using `rtacklayer::import` which returns a GRanges object 
with the coverage as "score" column. Bigwigs can be visualized in the IGV viewer as well in a memory-efficient fashion (unlike bedGraph).
If one aims to average bigwigs there are two options:
1) Use the script `AverageBigwig_v1.0.0.sh` in this repo.
2) Use `wiggletooms mean` which is cumbersome because it returns wig, which needs to bed converted back to bedGraph and then back to bigwig.
If one converts the wig directly to bigwig it will produce a very large files, much larger (in terms of Mb/Gb) than bigwig produced from bedGraph, no idea why.

<br>
<br>

#### AverageBigwig_v1.0.1.sh

This function takes two or more bigwigs and returns an averaged bigwig file.

```{bash}

------------------------------------------------------------------------------------------------------------------
AverageBigwig.sh
  
Usage: AverageBigwig.sh <inputfiles> <output> <chromsizes>
  
<inputfiles> : a space-delimited list of bigwigs, e.g. simply 
               obtained with <ls <whatever>*.bigwig>
               
<output>     : output file name
  
<chromsizes> : tab-delim text file indicating size of each chromosome
  
The recommended usage is:
$ ls *.bigwig | bash AverageBigwig.sh /dev/stdin output.bigwig chromsizes

------------------------------------------------------------------------------------------------------------------
```
  
One should parse the input files as a space-delimited string up front and then pipe this to the function as first argument.
Second argument is the output file and third one is a chromSizes file as the intermediate output is a bedGraph that needs to be converted
back to bigwig. Technically the script uses process substitution to feed the bigwigs into `bedtools unionbedg` like
`bedtools unionbedg <(bigwig1 /dev/stdout) <(bigwig2 /dev/stdout)` and then uses an `awk` one-liner to average this output which is 
basically the interval (chr-start-end) plus the coverage from each of the files. This intermediate result is then converted back to bigwig.
This is based on my [biostars post](https://www.biostars.org/p/329080/#329111) from a few years back and this again is extended based on this code snipped from
[Aaron Quinlan's Gist](https://gist.github.com/arq5x/5bdad2bd6d869ceca1ee).
The script does not have a notable memory footprint, one can easily average many bigwigs at once. 
It takes one thread for `bedtools` while both `mawk` (faster than `awk`) and the conversion tools don't need much CPU.

