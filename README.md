# NGS_Pipelines

Pipelines for preprocessing of ATAC/ChIP-seq, RNA-seq and scRNA-seq data. 

## Software

- the `environment.yml` contains a Linux conda environment (built on CentOS-7) with all required software

- a Docker image based on that environment is available from the [Docker Hub](https://hub.docker.com/r/atpoint/ngs_pipelines)

## Available Pipeline
  
Run any of the bash scripts without arguments or with `-h/--help` to see the help section with all available arguments.
  
### `RNAseq.sh`

**Deprecated:** Use this Nextflow pipeline instead => https://github.com/ATpoint/rnaseq_preprocess

### `scRNAseq.sh`
  
**Deprecated:** Use this Nextflow pipeline instead => https://github.com/ATpoint/sc_preprocess

### `DNAseq.sh`

**Deprecated:** Use this Nextflow pipeline instead => https://github.com/ATpoint/atac_chip_preprocess
 
### `Bam2Bigwig.sh`

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

### AverageBigwig.sh

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
$ ls *.bigwig | bash AverageBigwig.sh /dev/stdin output_averaged.bigwig chromsizes

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

