#!/bin/bash
#SBATCH --error=flagstat.%J.err
#SBATCH --output=flagstat.%J.out
#SBATCH --workdir=/data/htp/A07/RRBS/AML_triplets/bsmap
#SBATCH --mem=40000

for i in /data/htp/A07/RRBS/AML_triplets/bsmap/*_bsmap_primary_sorted.bam
do
echo $i >> flagstat_output.txt
samtools flagstat $i  >> flagstat_output.txt

echo $i >> number_mapped_reads.txt
samtools view -F 0x04 -c $i  >> number_mapped_reads.txt
done
