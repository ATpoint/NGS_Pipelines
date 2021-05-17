# NGS_Pipelines

Bash-based pipelines for processing of ATAC/ChIP-seq, RNA-seq and scRNA-seq data.
Nextflow replacements will follow soon.


## Software Installation

A Docker container with all necessary software is available from the [Docker Hub](https://hub.docker.com/repository/docker/atpoint/phd_project) based on this [Github repo](https://github.com/ATpoint/phd_project_docker). Inside the container one has to run the following commands to activate the conda environment. On the HPC one can pull the image and convert to the Singularity `sif` format using as shown below. The environment export is shown [here](#conda-environment-export).

```

#/ activate env:
eval "$(conda shell.bash hook)"
conda activate Pipelines

#/ on HPC pull and convert to SIF:
singularity pull docker://atpoint90/phd_project:1.0.0

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

## Conda environment export

```

$ conda env export --name Pipelines

name: Pipelines
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=1_gnu
  - _r-mutex=1.0.1=anacondar_1
  - alsa-lib=1.2.3=h516909a_0
  - bedtools=2.29.2=hc088bd4_0
  - binutils_impl_linux-64=2.35=h18a2f87_9
  - binutils_linux-64=2.35=hc3fd857_29
  - bioconductor-edger=3.32.0=r40h5f743cb_0
  - bioconductor-limma=3.46.0=r40h037d062_0
  - bowtie2=2.4.2=py38h1c8e9b9_1
  - brotlipy=0.7.0=py38h8df0ef7_1001
  - bwidget=1.9.14=0
  - bzip2=1.0.8=h516909a_3
  - c-ares=1.16.1=h516909a_3
  - ca-certificates=2021.4.13=h06a4308_1
  - cairo=1.16.0=hcf35c78_1003
  - certifi=2020.12.5=py38h578d9bd_1
  - cffi=1.14.3=py38h1bdcb99_1
  - chardet=3.0.4=py38h924ce5b_1008
  - click=7.1.2=pyh9f0ad1d_0
  - cmake=3.18.4=h1f3970d_0
  - coloredlogs=14.0=py38h32f6830_2
  - colormath=3.0.0=py_2
  - coreutils=8.31=h516909a_0
  - cryptography=3.2=py38hb23e4d4_0
  - curl=7.71.1=he644dc0_8
  - cutadapt=2.10=py38h0213d0e_1
  - cycler=0.10.0=py_2
  - decorator=4.4.2=py_0
  - dnaio=0.4.2=py38h0213d0e_1
  - expat=2.2.9=he1b5a44_2
  - fastqc=0.11.9=0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - fontconfig=2.13.1=h86ecdb6_1001
  - freetype=2.10.4=he06d7ca_0
  - fribidi=1.0.10=h516909a_0
  - future=0.18.2=py38h32f6830_2
  - gcc_impl_linux-64=7.5.0=hda68d29_13
  - gcc_linux-64=7.5.0=he2a3fca_29
  - genrich=0.6.1=h5bf99c6_1
  - gettext=0.19.8.1=hf34092f_1004
  - gfortran_impl_linux-64=7.5.0=hfca37b7_17
  - gfortran_linux-64=7.5.0=ha081f1e_29
  - ghostscript=9.53.3=h58526e2_2
  - giflib=5.2.1=h516909a_2
  - glib=2.66.2=he1b5a44_0
  - graphite2=1.3.13=he1b5a44_1001
  - gsl=2.6=h294904e_0
  - gxx_impl_linux-64=7.5.0=h64c220c_13
  - gxx_linux-64=7.5.0=h547f3ba_29
  - harfbuzz=2.4.0=h9f30f68_3
  - htslib=1.11=hd3b49d5_0
  - humanfriendly=8.2=py38h32f6830_1
  - icu=64.2=he1b5a44_1
  - idna=2.10=pyh9f0ad1d_0
  - importlib-metadata=2.0.0=py_1
  - jemalloc=5.2.1=he1b5a44_4
  - jinja2=2.11.2=pyh9f0ad1d_0
  - jpeg=9d=h516909a_0
  - kernel-headers_linux-64=2.6.32=h77966d4_13
  - kiwisolver=1.2.0=py38hbf85e49_1
  - krb5=1.17.1=hfafb76e_3
  - lcms2=2.11=hbd6801e_0
  - ld_impl_linux-64=2.35=h769bd43_9
  - libblas=3.9.0=2_openblas
  - libcblas=3.9.0=2_openblas
  - libcurl=7.71.1=hcdd3856_8
  - libdeflate=1.6=h516909a_0
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=h516909a_1
  - libffi=3.2.1=he1b5a44_1007
  - libgcc=7.2.0=h69d50b8_2
  - libgcc-devel_linux-64=9.3.0=hfd08b2a_17
  - libgcc-ng=9.3.0=h5dbcf3e_17
  - libgfortran-ng=7.5.0=hae1eefd_17
  - libgfortran4=7.5.0=hae1eefd_17
  - libglib=2.66.2=h0dae87d_0
  - libgomp=9.3.0=h5dbcf3e_17
  - libiconv=1.16=h516909a_0
  - liblapack=3.9.0=2_openblas
  - libnghttp2=1.41.0=h8cfc5f6_2
  - libopenblas=0.3.12=pthreads_hb3c22a3_1
  - libpng=1.6.37=hed695b0_2
  - libssh2=1.9.0=hab1572f_5
  - libstdcxx-devel_linux-64=9.3.0=h4084dd6_17
  - libstdcxx-ng=9.3.0=h2ae2ef3_17
  - libtiff=4.1.0=hc7e4089_6
  - libuuid=2.32.1=h14c3975_1000
  - libuv=1.40.0=hd18ef5c_0
  - libwebp-base=1.1.0=h516909a_3
  - libxcb=1.13=h14c3975_1002
  - libxml2=2.9.10=hee79883_0
  - libxslt=1.1.33=h31b3aaa_0
  - lz4-c=1.9.2=he1b5a44_3
  - lzstring=1.0.4=py_1001
  - macs2=2.2.7.1=py38h0213d0e_1
  - make=4.3=hd18ef5c_1
  - markdown=3.3.3=pyh9f0ad1d_0
  - markupsafe=1.1.1=py38h8df0ef7_2
  - matplotlib-base=3.3.2=py38h4d1ce4f_1
  - mawk=1.3.4=h516909a_3
  - meme=5.3.0=py38pl526hc1f1133_0
  - mpi=1.0=openmpi
  - multiqc=1.9=py_1
  - mysql-connector-c=6.1.11=h6eb9d5d_1007
  - ncurses=6.2=he1b5a44_2
  - networkx=2.5=py_0
  - numpy=1.19.2=py38hf89b668_1
  - olefile=0.46=pyh9f0ad1d_1
  - openjdk=11.0.8=hacce0ff_0
  - openmpi=4.0.5=ha4a8674_4
  - openssl=1.1.1k=h7f98852_0
  - pango=1.42.4=h7062337_4
  - parallel=20200922=0
  - pcre=8.44=he1b5a44_0
  - pcre2=10.35=h279444b_1
  - perl=5.26.2=h36c2ea0_1008
  - perl-app-cpanminus=1.7044=pl526_1
  - perl-base=2.23=pl526_1
  - perl-carp=1.38=pl526_3
  - perl-cgi=4.44=pl526h14c3975_1
  - perl-common-sense=3.74=pl526_2
  - perl-constant=1.33=pl526_1
  - perl-dbi=1.642=pl526_0
  - perl-encode=2.88=pl526_1
  - perl-exporter=5.72=pl526_1
  - perl-extutils-makemaker=7.36=pl526_1
  - perl-file-path=2.16=pl526_0
  - perl-file-temp=0.2304=pl526_2
  - perl-file-which=1.23=pl526_0
  - perl-html-parser=3.72=pl526h6bb024c_5
  - perl-html-tagset=3.20=pl526_3
  - perl-html-template=2.97=pl526_1
  - perl-html-tree=5.07=pl526_1
  - perl-json=4.02=pl526_0
  - perl-json-xs=2.34=pl526h6bb024c_3
  - perl-log-log4perl=1.49=pl526_0
  - perl-math-cdf=0.1=pl526h14c3975_5
  - perl-parent=0.236=pl526_1
  - perl-scalar-list-utils=1.52=pl526h516909a_0
  - perl-types-serialiser=1.0=pl526_2
  - perl-xml-namespacesupport=1.12=pl526_0
  - perl-xml-parser=2.44_01=pl526ha1d75be_1002
  - perl-xml-sax=1.02=pl526_0
  - perl-xml-sax-base=1.09=pl526_0
  - perl-xml-sax-expat=0.51=pl526_3
  - perl-xml-simple=2.25=pl526_1
  - perl-xsloader=0.24=pl526_0
  - perl-yaml=1.29=pl526_0
  - picard=2.23.8=0
  - pigz=2.3.4=hed695b0_1
  - pillow=8.0.1=py38h9776b28_0
  - pip=20.2.4=py_0
  - pixman=0.38.0=h516909a_1003
  - pthread-stubs=0.4=h14c3975_1001
  - pycparser=2.20=pyh9f0ad1d_2
  - pyopenssl=19.1.0=py_1
  - pyparsing=2.4.7=pyh9f0ad1d_0
  - pysocks=1.7.1=py38h924ce5b_2
  - python=3.8.6=h852b56e_0_cpython
  - python-dateutil=2.8.1=py_0
  - python_abi=3.8=1_cp38
  - pyyaml=5.3.1=py38h8df0ef7_1
  - r-base=4.0.1=h95c6c4b_0
  - r-lattice=0.20_41=r40hcdcec82_2
  - r-locfit=1.5_9.4=r40hcdcec82_1
  - r-rcpp=1.0.4.6=r40h0357c0b_1
  - readline=8.0=he28a2e2_2
  - requests=2.24.0=pyh9f0ad1d_0
  - rhash=1.3.6=h516909a_1001
  - salmon=1.3.0=hf69c8f4_0
  - samblaster=0.1.26=hc9558a2_0
  - samtools=1.11=h6270b1f_0
  - sed=4.8=hbfbb72e_0
  - seqtk=1.3=hed695b0_2
  - setuptools=49.6.0=py38h924ce5b_2
  - simplejson=3.17.2=py38h1e0a361_1
  - six=1.15.0=pyh9f0ad1d_0
  - spectra=0.0.11=py_1
  - sqlite=3.33.0=h4cf870e_1
  - subread=2.0.1=hed695b0_0
  - sysroot_linux-64=2.12=h77966d4_13
  - tabix=0.2.6=ha92aebf_0
  - tbb=2020.2=hc9558a2_0
  - tk=8.6.10=hed695b0_1
  - tktable=2.10=h555a92e_3
  - tornado=6.0.4=py38h1e0a361_2
  - ucsc-bedgraphtobigwig=377=h446ed27_1
  - ucsc-bigwigtobedgraph=377=h446ed27_4
  - urllib3=1.25.11=py_0
  - wheel=0.35.1=pyh9f0ad1d_0
  - xopen=0.9.0=py38h32f6830_1
  - xorg-fixesproto=5.0=h14c3975_1002
  - xorg-inputproto=2.3.2=h14c3975_1002
  - xorg-kbproto=1.0.7=h14c3975_1002
  - xorg-libice=1.0.10=h516909a_0
  - xorg-libsm=1.2.3=h84519dc_1000
  - xorg-libx11=1.6.12=h516909a_0
  - xorg-libxau=1.0.9=h14c3975_0
  - xorg-libxdmcp=1.1.3=h516909a_0
  - xorg-libxext=1.3.4=h516909a_0
  - xorg-libxfixes=5.0.3=h516909a_1004
  - xorg-libxi=1.7.10=h516909a_0
  - xorg-libxrender=0.9.10=h516909a_1002
  - xorg-libxtst=1.2.3=h516909a_1002
  - xorg-recordproto=1.14.2=h516909a_1002
  - xorg-renderproto=0.11.1=h14c3975_1002
  - xorg-xextproto=7.3.0=h14c3975_1002
  - xorg-xproto=7.0.31=h14c3975_1007
  - xz=5.2.5=h516909a_1
  - yaml=0.2.5=h516909a_0
  - zipp=3.4.0=py_0
  - zlib=1.2.11=h516909a_1010
  - zstd=1.4.5=h6597ccf_2

```
