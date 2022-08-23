#!/bin/bash


slurms="/data/htp/A07/RRBS/AML_triplets/TestData/bsmap/slurm"

for prefix in `cut -f2 /data/htp/A07/RRBS/AML_triplets/TestData/deML/index.txt`
do
# MAKING THE HEADER
echo '#!/bin/bash' > $slurms/$prefix.bsmap.sh
echo '#SBATCH --workdir='$slurms >> $slurms/$prefix.bsmap.sh
echo '#SBATCH --error='$prefix'.%J.err' >>$slurms/$prefix.bsmap.sh
echo '#SBATCH --output='$prefix'.%J.out' >> $slurms/$prefix.bsmap.sh
echo '#SBATCH --cpus-per-task=3' >> $slurms/$prefix.bsmap.sh
echo '#SBATCH --mem=1000' >> $slurms/$prefix.bsmap.sh

# mapping
echo "bsmap -a /data/htp/A07/RRBS/cutadapt/test_"$prefix"_cutadapt.r1.fq.gz \
 -b /data/htp/A07/RRBS/AML_triplets/cutadapt/test_"$prefix"_cutadapt.r2.fq.gz \
 -d /data/ngs/genomes/Human/hg38/hg38.fa -o /data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix".bsmap_output.bam \
 -D C-CGG -u -p 10" >> $slurms/$prefix.bsmap.sh

# primary alignment
echo "samtools view -F 0x100 -b /data/htp/A07/RRBS/AML_triplets/bsmap/$prefix.bsmap_output.bam >> \
/data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary.bam" >> $slurms/$prefix.bsmap.sh

# sort bam file
echo "samtools sort /data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary.bam > \
/data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary_sorted.bam" >> $slurms/$prefix.bsmap.sh

# mark duplicates
echo "picard MarkDuplicates \
     I= /data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_primary_sorted.bam \
     O=/data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_markdup.bam \
     M=/data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_marked_dup_metrics.txt" >> $slurms/$prefix.bsmap.sh

# filter out x & y chromosome
echo "samtools idxstats /data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary_sorted.bam | cut -f 1 | grep -v -E 'chrX|chrY'  | \
xargs samtools view -b /data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary_sorted.bam > \
/data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary_sorted_filter.bam" >> $slurms/$prefix.bsmap.sh

# index bam
echo "samtools index /data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_primary_sorted_filter.bam" >> $slurms/$prefix.bsmap.sh
echo "samtools index /data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_markdup.bam" >> $slurms/$prefix.bsmap.sh

# make methratio file
echo "python /data/htp/A07/RRBS/common/methratio.py /data/htp/A07/RRBS/AML_triplets/bsmap/"$prefix"_bsmap_primary_sorted_filter.bam \
-o /data/htp/A07/RRBS/AML_triplets/methratio/filtered/"$prefix".methratio.txt -d /data/ngs/genomes/Human/hg38/hg38.fa"  >> $slurms/$prefix.bsmap.sh


sbatch $slurms/$prefix.bsmap.sh

done
