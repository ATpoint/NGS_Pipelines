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

**Deprecated:** Use this https://github.com/ATpoint/bam2bedgraph, still bash-based but a more cleaned-up version that
is simpler and includes th averaging script from below.

### `AverageBigwig.sh

**Deprecated**: See https://github.com/ATpoint/bam2bedgraph.

