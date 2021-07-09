#!/bin/bash

#/ Sort files into folders:

function exist_move {

  if [[ $(find . -maxdepth 1 -name "${1}" | wc -l) > 0 ]]; then
    mkdir -p ${2}
    mv $1 ./${2}/
  fi
}

rm *.fastq.gz

exist_move "*raw.bam*" "bam_raw"; exist_move "*raw.flag*" "*raw.flagstat"

exist_move "*_dedup.bam*" "./bam_dedup"; exist_move "*_dedup.flag*" "./bam_dedup"

exist_move "*.bigwig" "./bigwig"

exist_move "*fastqc.*" "./fastqc"

exist_move "*peaks*" "./macs2"

exist_move "*.log" "./logs"

exist_move "*InsertSize*" "./insert_sizes"
