#!/bin/bash
#SBATCH --error=qc.%J.err
#SBATCH --output=qc.%J.out
#SBATCH --workdir=/data/htp/A07/RRBS/AML_triplets/TotalData/QC/
#SBATCH --cpus-per-task=4
#SBATCH --mem=3000



srun  Rscript /home/weser/AML_triplets/RRBS/TotalData/QC.R
