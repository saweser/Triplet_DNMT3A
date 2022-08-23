#!/bin/bash
#SBATCH --error=deML.%J.err
#SBATCH --output=deML.%J.out
#SBATCH --workdir=/data/htp/A07/RRBS/deML
#SBATCH --mem=100000





deML -i /data/htp/A07/RRBS/deML/index_491.txt -f /data/htp/A07/RRBS/fastq_files/test_lib/test_lib_R1.fq.gz -r /data/htp/A07/RRBS/fastq_files/test_lib/test_lib_R4.fq.gz -if1 /data/htp/A07/RRBS/fastq_files/test_lib/test_lib_R2.fq.gz -o /data/htp/A07/RRBS/deML/test -s summary.txt

# you should have more than 80% of reads successfully demultiplexed