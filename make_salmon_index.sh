#!/bin/bash

#/ Create indices for human and mouse using GENCODE annotations:

#------------------------------------------------------------------------------------------------------------------

###/ Human:

#/ Download from GENCODE, remove the whitespace in the header, compress with bgzip:
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/GRCh38.primary_assembly.genome.fa.gz && \
zcat GRCh38.primary_assembly.genome.fa.gz \
| mawk '$1 ~ /^>/ {split($1,a," "); print a[1]; next} { print }' \
| bgzip -@ 8 > foo && \
mv foo GRCh38.primary_assembly.genome.fa.gz
 
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.transcripts.fa.gz && \
zcat gencode.v35.transcripts.fa.gz | bgzip -@ 8 > foo && \
mv foo gencode.v35.transcripts.fa.gz

ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz && \
zcat gencode.v35.annotation.gtf.gz | bgzip -@ 8 > foo && \
mv foo gencode.v35.annotation.gtf.gz

#/ Indexing:
mkdir -p ./Index/index_salmon/

#/ version:
salmon --version

#/ essential part of index creation:
salmon swim

#/ get chromosome names for decoying:
zgrep '^>' GRCh38.primary_assembly.genome.fa.gz \
| cut -d ">" -f 2 > ./Index/index_salmon/decoys.txt

cat gencode.v35.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > ./Index/index_salmon/gentrome.fa.gz 

#/ build default index:
salmon index  \
--transcripts ./Index/index_salmon/gentrome.fa.gz \
--threads 16 --gencode --kmerLen 31 --sparse \
--decoys ./Index/index_salmon/decoys.txt --index ./Index/index_salmon/kmer31
