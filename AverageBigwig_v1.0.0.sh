#!/bin/bash

#/ A simple wrapper to average bigwig files via bedtools unionbedg and awk.

usage(){
  echo '
  
  --------------------------------------------------------------------------------------------------------
  
  AverageBigwig.sh
  
  Usage: AverageBigwig <inputfiles> <output> <chromsizes>
  
  <inputfiles> : a space-delimited list of bigwigs, e.g. simply 
                 obtained with <ls <whatever>*.bigwig>
               
  <output>     : output file name
  
  <chromsizes> : tab-delim text file indicating size of each chromosome
  
  The recommended usage is:
  $ ls *.bigwig | bash AverageBigwig.sh /dev/stdin output.bigwig chromsizes
  
  --------------------------------------------------------------------------------------------------------
               
  '
}; if [[ -z "$1" ]] || [[ $1 == -h ]] || [[ $1 == --help ]]; then usage; exit; fi

#---------------------------------------------------------------------------------------------------------

if [ $# -lt 3 ]; then
    echo '[ERROR] Command must have three arguments!' && usage
    exit 1
fi

#---------------------------------------------------------------------------------------------------------


if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

function PathCheck {
  
  if [[ $(command -v $1 | wc -l | xargs) == 0 ]]; then 
  echo ${1} >> missing_tools.txt
  fi
  
}; export -f PathCheck

Tools=(bedGraphToBigWig bedtools bigWigToBedGraph mawk)

#/ if that argument is set check for tools, then exit:
if [[ -e missing_tools.txt ]]; then rm missing_tools.txt; fi

for i in $(echo ${Tools[*]}); do PathCheck $i; done

if [[ -e missing_tools.txt ]]; then 
  echo '[ERROR] Tools missing in PATH, see missing_tools.txt'
  exit 1
fi 

#---------------------------------------------------------------------------------------------------------

#/ Function:
Input=$(awk '{print $0}' "${1}")
Output="${2}"
ChromSizes="${3}"
  
#/ Parse input files into a process-substitution string:
Bigwigs=$(echo ${Input} \
          | tr " " "\n" \
          | while read p; do echo "<(bigWigToBedGraph ${p} /dev/stdout)"; done \
          < /dev/stdin)
          
#/ Put it together:            
eval bedtools unionbedg -i $Bigwigs \
| mawk 'OFS="\t" {sum=0; for (col=4; col<=NF; col++) sum += $col; print $1, $2, $3, sum/(NF-4+1); }' \
> "${Output}".tmp.bedGraph

if [[ -e "${Output}".tmp.bedGraph ]]; then 

  bedGraphToBigWig "${Output}".tmp.bedGraph "${ChromSizes}" "${Output}" && rm "${Output}".tmp.bedGraph \
  || exit 1
  
fi  
  
#---------------------------------------------------------------------------------------------------------
