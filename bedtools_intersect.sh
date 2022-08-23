#!/bin/bash
#SBATCH --error=methylkit.%J.err
#SBATCH --output=methylkit.%J.out
#SBATCH --workdir=/data/htp/A07/RRBS
#SBATCH --mem=4000

bedtools intersect -b /data/htp/A07/RRBS/tables/downsampled_all_Cs_18x.bed -a /data/htp/A07/RRBS/tables/anno_CpGs_hsap_ucsc.bed -header > /data/htp/A07/RRBS/tables/all_islands_downsampled_all_Cs_18x.bed &
bedtools intersect -b /data/htp/A07/RRBS/tables/downsampled_all_Cs_15x.bed -a /data/htp/A07/RRBS/tables/anno_CpGs_hsap_ucsc.bed -header > /data/htp/A07/RRBS/tables/all_islands_downsampled_all_Cs_15x.bed &
bedtools intersect -b /data/htp/A07/RRBS/tables/downsampled_all_Cs_10x.bed -a /data/htp/A07/RRBS/tables/anno_CpGs_hsap_ucsc.bed -header > /data/htp/A07/RRBS/tables/all_islands_downsampled_all_Cs_10x.bed &
bedtools intersect -b /data/htp/A07/RRBS/tables/downsampled_all_Cs_5x.bed -a /data/htp/A07/RRBS/tables/anno_CpGs_hsap_ucsc.bed -header > /data/htp/A07/RRBS/tables/all_islands_downsampled_all_Cs_5x.bed &
bedtools intersect -b /data/htp/A07/RRBS/tables/downsampled_all_Cs_2x.bed -a /data/htp/A07/RRBS/tables/anno_CpGs_hsap_ucsc.bed -header > /data/htp/A07/RRBS/tables/all_islands_downsampled_all_Cs_2x.bed &
bedtools intersect -b /data/htp/A07/RRBS/tables/downsampled_all_Cs_1x.bed -a /data/htp/A07/RRBS/tables/anno_CpGs_hsap_ucsc.bed -header > /data/htp/A07/RRBS/tables/all_islands_downsampled_all_Cs_1x.bed &

