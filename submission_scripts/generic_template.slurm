#!/bin/bash

#SBATCH --export=NONE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --mem=80G
#SBATCH --partition=normal
#SBATCH --time=24:00:00
#SBATCH --job-name=jobname
#SBATCH --output=jobname.log
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
  
#/ Unstranded RNA-seq in mice:
./RNAseq_v1.0.2.sh \
  --idx /scratch/tmp/a_toen03/Genomes/GRCm38/Index/index_salmon/kmer31/ \
  --mode paired \
  --threads 10 \
  --njobs 3 \
  --trim \
  --trimthreads 6 \
  --trimjobs 4 \
  --libtype IU
  
#/ scRNA-seq 10X-V3 in mice:
Idxdir="/scratch/tmp/a_toen03/Genomes/GRCm38/Index/index_alevin/"
./scRNAseq_v1.0.0.sh \
  --idx ${Idxdir}/gencode.vM25.annotation.expanded.sidx/ \
  --tgmap ${Idxdir}/gencode.vM25.annotation.expanded.tx2gene.tsv \
  --mtgenes ${Idxdir}/mtRNA.txt \
  --rrnagenes ${Idxdir}/rRNA.txt \
  --threads 8 \
  --njobs 4 \
  --additional "--numCellBootstraps 100 --dumpFeatures" 
    
#/ Run fastqc on everything:  
ls *.fastq.gz | parallel -j 30 "fastqc {}"

#/ Summary report:
multiqc .
