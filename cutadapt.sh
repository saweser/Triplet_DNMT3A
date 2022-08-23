#!/bin/bash
#SBATCH --error=cutadapt.%J.err
#SBATCH --output=cutadapt.%J.out
#SBATCH --workdir=/data/htp/A07/RRBS/cutadapt
#SBATCH --mem=100000

for prefix in `cut -f2 /data/htp/A07/RRBS/deML/index.txt`
do
cutadapt -a NNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
-A NNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --quality-cutoff=20 \
--minimum-length=35 --overlap=1 -o /data/htp/A07/RRBS/cutadapt/test_"$prefix"_cutadapt.r1.fq.gz \
-p /data/htp/A07/RRBS/cutadapt/test_"$prefix"_cutadapt.r2.fq.gz \
/data/htp/A07/RRBS/deML/test_"$prefix"_r1.fq.gz \
/data/htp/A07/RRBS/deML/test_"$prefix"_r2.fq.gz
done
