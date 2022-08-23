#!/bin/bash
#SBATCH --workdir=/data/htp/A07/RRBS/bsmap/methratio/filtered_secondary_alignment

slurms="/data/htp/A07/RRBS/bsmap/methratio/filtered_secondary_alignment/CpGs/slurm"


for file in `ls *.txt | sed 's/.txt//g'`;
do 

# MAKING THE HEADER
echo '#!/bin/bash' > $slurms/$file.methratio.sh
echo '#SBATCH --workdir='$slurms >> $slurms/$file.methratio.sh
echo '#SBATCH --error='$file'.%J.err' >>$slurms/$file.methratio.sh
echo '#SBATCH --output='$file'.%J.out' >> $slurms/$file.methratio.sh
echo '#SBATCH --cpus-per-task=10' >> $slurms/$file.methratio.sh
echo '#SBATCH --mem=3000' >> $slurms/$file.methratio.sh

awk '{if($4=="CG")print;}' /data/htp/A07/RRBS/bsmap/methratio/filtered_secondary_alignment/$file.txt >> /data/htp/A07/RRBS/bsmap/methratio/filtered_secondary_alignment/CpGs/CpGs_$file.txt
sbatch $slurms/$file.methratio.sh

done
