#!/bin/bash

#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --mem=80G
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --job-name=DNAseq
#SBATCH --output=DNAseq_v1.0.x.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=foo@bar.de

eval "$(conda shell.bash hook)"
conda activate Pipelines

ulimit -u 50000

#/ Example for a paired-end ATAC-seq experiment in mice:
./DNAseq_v1.0.0.sh \
  --genome /scratch/tmp/a_toen03/Genomes/GRCm38/Index/index_bowtie2/GRCm38.primary_assembly.genome.fa \
  --jobnumber 4 \
  --threads 14 \
  --format fq_pe \
  --atacseq

#/ Single-end ChIP-seq data:
./DNAseq_v1.0.0.sh \
  --genome /scratch/tmp/a_toen03/Genomes/GRCm38/Index/index_bowtie2/GRCm38.primary_assembly.genome.fa \
  --jobnumber 4 \
  --threads 16 \
  --format fq_se
  
#/ Run fastqc on everything:  
ls *.fastq.gz | parallel -j 30 "fastqc {}"

#/ Summary report:
multiqc .
