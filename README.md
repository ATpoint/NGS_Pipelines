# NGS_Pipelines

Bash-based pipelines for processing of ATAC/ChIP-seq, RNA-seq and scRNA-seq data.
Not intended for public use, therefore use at your own risk, without any warranty.

## Software Installation
It is recommended to install software via the `miniconda` package manager into a separate environment.
A Linux system and `miniconda3` is assumed, on Mac it should work as well, but not extensively tested.
After installing `conda` one should change `.condarc` to:

```{bash}

channels:
  - conda-forge
  - bioconda
  - defaults

```

...and then create the environment/install the software.

```

#/ Install conda if not already done:
if [[ "$(uname)" == "Darwin" ]]; then
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
elif [[ "$(uname)" == "Linux" ]]; then
  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
fi

#/ If on Mac and brew is the main pkg manager then run this to avoid auto-activate of conda base:
echo 'auto_activate_base: false' >> ~/.condarc  

#/ Create environment and install necessary tools:
conda create --name Pipelines
conda activate Pipelines

Tools=(bedtools=2.29.2 bioconductor-edger=3.32.0 bowtie2=2.4.2 coreutils cutadapt \
       fastqc subread macs2=2.2.7.1 mawk multiqc r-base=4.0.1 parallel picard pigz \
       salmon=1.3.0 samblaster samtools=1.11 seqtk tabix ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph)
        
printf '%s\n' "${Tools[@]}" > install_software.txt

conda install --file install_software.txt

```

In order to activate the envir from inside a SLURM script add these lines to the top of the script as 
suggested at [this Github issue](https://github.com/conda/conda/issues/7980).

```
eval "$(conda shell.bash hook)"
conda activate Pipelines
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

#### AverageBigwig_v1.0.0.sh

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

#### Software versions

This is the output of `conda list | grep -v '#' | awk '{print $1"="$2}'` which are the versions of software to run the pipelines on on our HPC. Saving this to a file `software_versioned.txt` should enable to run `conda install --file software_versioned.txt` to reproduce our environment.

```{bash}

_libgcc_mutex=0.1
_openmp_mutex=4.5
_r-mutex=1.0.1
alsa-lib=1.2.3
bedtools=2.29.2
binutils_impl_linux-64=2.35
binutils_linux-64=2.35
bioconductor-edger=3.32.0
bioconductor-limma=3.46.0
bowtie2=2.4.2
brotlipy=0.7.0
bwidget=1.9.14
bzip2=1.0.8
c-ares=1.16.1
ca-certificates=2020.6.20
cairo=1.16.0
certifi=2020.6.20
cffi=1.14.3
chardet=3.0.4
click=7.1.2
cmake=3.18.4
coloredlogs=14.0
colormath=3.0.0
coreutils=8.31
cryptography=3.2
curl=7.71.1
cutadapt=2.10
cycler=0.10.0
decorator=4.4.2
dnaio=0.4.2
expat=2.2.9
fastqc=0.11.9
font-ttf-dejavu-sans-mono=2.37
fontconfig=2.13.1
freetype=2.10.4
fribidi=1.0.10
future=0.18.2
gcc_impl_linux-64=7.5.0
gcc_linux-64=7.5.0
gettext=0.19.8.1
gfortran_impl_linux-64=7.5.0
gfortran_linux-64=7.5.0
giflib=5.2.1
glib=2.66.2
graphite2=1.3.13
gsl=2.6
gxx_impl_linux-64=7.5.0
gxx_linux-64=7.5.0
harfbuzz=2.4.0
htslib=1.11
humanfriendly=8.2
icu=64.2
idna=2.10
importlib-metadata=2.0.0
jemalloc=5.2.1
jinja2=2.11.2
jpeg=9d
kernel-headers_linux-64=2.6.32
kiwisolver=1.2.0
krb5=1.17.1
lcms2=2.11
ld_impl_linux-64=2.35
libblas=3.9.0
libcblas=3.9.0
libcurl=7.71.1
libdeflate=1.6
libedit=3.1.20191231
libev=4.33
libffi=3.2.1
libgcc=7.2.0
libgcc-devel_linux-64=9.3.0
libgcc-ng=9.3.0
libgfortran-ng=7.5.0
libgfortran4=7.5.0
libglib=2.66.2
libgomp=9.3.0
libiconv=1.16
liblapack=3.9.0
libnghttp2=1.41.0
libopenblas=0.3.12
libpng=1.6.37
libssh2=1.9.0
libstdcxx-devel_linux-64=9.3.0
libstdcxx-ng=9.3.0
libtiff=4.1.0
libuuid=2.32.1
libuv=1.40.0
libwebp-base=1.1.0
libxcb=1.13
libxml2=2.9.10
lz4-c=1.9.2
lzstring=1.0.4
macs2=2.2.7.1
make=4.3
markdown=3.3.3
markupsafe=1.1.1
matplotlib-base=3.3.2
mawk=1.3.4
multiqc=1.9
mysql-connector-c=6.1.11
ncurses=6.2
networkx=2.5
numpy=1.19.2
olefile=0.46
openjdk=11.0.8
openssl=1.1.1h
pango=1.42.4
parallel=20200922
pcre=8.44
pcre2=10.35
perl=5.30.3
picard=2.23.8
pigz=2.3.4
pillow=8.0.1
pip=20.2.4
pixman=0.38.0
pthread-stubs=0.4
pycparser=2.20
pyopenssl=19.1.0
pyparsing=2.4.7
pysocks=1.7.1
python=3.8.6
python-dateutil=2.8.1
python_abi=3.8
pyyaml=5.3.1
r-base=4.0.1
r-lattice=0.20_41
r-locfit=1.5_9.4
r-rcpp=1.0.4.6
readline=8.0
requests=2.24.0
rhash=1.3.6
salmon=1.3.0
samblaster=0.1.26
samtools=1.11
sed=4.8
seqtk=1.3
setuptools=49.6.0
simplejson=3.17.2
six=1.15.0
spectra=0.0.11
sqlite=3.33.0
subread=2.0.1
sysroot_linux-64=2.12
tabix=0.2.6
tbb=2020.2
tk=8.6.10
tktable=2.10
tornado=6.0.4
ucsc-bedgraphtobigwig=377
ucsc-bigwigtobedgraph=377
urllib3=1.25.11
wheel=0.35.1
xopen=0.9.0
xorg-fixesproto=5.0
xorg-inputproto=2.3.2
xorg-kbproto=1.0.7
xorg-libice=1.0.10
xorg-libsm=1.2.3
xorg-libx11=1.6.12
xorg-libxau=1.0.9
xorg-libxdmcp=1.1.3
xorg-libxext=1.3.4
xorg-libxfixes=5.0.3
xorg-libxi=1.7.10
xorg-libxrender=0.9.10
xorg-libxtst=1.2.3
xorg-recordproto=1.14.2
xorg-renderproto=0.11.1
xorg-xextproto=7.3.0
xorg-xproto=7.0.31
xz=5.2.5
yaml=0.2.5
zipp=3.4.0
zlib=1.2.11
zstd=1.4.5

```
<chromsizes> : tab-delim text file indicating size of each chromosome
