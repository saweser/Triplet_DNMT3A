#!/bin/bash
#SBATCH --error=fastqc.%J.err
#SBATCH --output=fastqc.%J.out
#SBATCH --workdir=/data/htp/A07/RRBS
#SBATCH --mem=100000



srun  Rscript /data/htp/A07/RRBS/fastqc.R


