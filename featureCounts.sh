#!/bin/bash
#SBATCH --error=featureCounts.%J.err
#SBATCH --output=featureCounts.%J.out
#SBATCH --workdir=/data/htp/A07/AML_triplets/expression
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=6



srun  Rscript /data/htp/A07/AML_triplets/feature_counts.R
