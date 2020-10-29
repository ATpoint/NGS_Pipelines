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
#SBATCH --mail-user=a_toen03@uni-muenster.de

eval "$(conda shell.bash hook)"
conda activate Pipelines
