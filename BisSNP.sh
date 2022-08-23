#!/bin/bash

slurms="/data/htp/A07/RRBS/AML_triplets/TestData/BisSNP/slurm"
n=0
for prefix in `cut -f2 /data/htp/A07/RRBS/AML_triplets/TestData/deML/index.txt`
do
n=$( expr $n + 1)
# MAKING THE HEADER
echo '#!/bin/bash' > $slurms/$prefix.bissnp.sh
echo '#SBATCH --workdir='$slurms >> $slurms/$prefix.bissnp.sh
echo '#SBATCH --error='$prefix'.%J.err' >>$slurms/$prefix.bissnp.sh
echo '#SBATCH --output='$prefix'.%J.out' >> $slurms/$prefix.bissnp.sh
echo '#SBATCH --cpus-per-task=3' >> $slurms/$prefix.bissnp.sh
echo '#SBATCH --mem=1000' >> $slurms/$prefix.bissnp.sh


# add read groups (GATK needs those)
echo "picard AddOrReplaceReadGroups \
I=/data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_primary_sorted_filter.bam \
O=/data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_primary_sorted_filter_rg.bam \
RGID=$n \
RGLB=lib1 \
RGPL=illumina \
RGPU=run \
RGSM=$prefix \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
SORT_ORDER=coordinate" >> $slurms/$prefix.bissnp.sh


echo "java -Xmx4g -jar /home/weser/BisSNP/BisSNP-1.0.0.jar \
-R /data/ngs/genomes/Human/hg38/hg38.fa \
-T BisulfiteGenotyper \
-I /data/htp/A07/RRBS/AML_triplets/TestData/bsmap/"$prefix"_bsmap_primary_sorted_filter_rg.bam \
-D /data/ngs/genomes/Human/hg38/gatk_bundle/hg38bundle/dbsnp_144.hg38.chr.vcf.gz \
-vfn1 /data/htp/A07/RRBS/AML_triplets/TestData/BisSNP/"$prefix"_cpg.raw.vcf \
-vfn2 /data/htp/A07/RRBS/AML_triplets/TestData/BisSNP/"$prefix"_snp.raw.vcf \
-L /data/ngs/genomes/Human/hg38/gatk_bundle/hg38bundle/wholegenome_interval_list.bed" >> $slurms/$prefix.bissnp.sh

# convert to bed format with perl script and filter only CpGs with CG
echo "perl /home/weser/BisSNP/vcf2bed.pl /data/htp/A07/RRBS/AML_triplets/TestData/BisSNP/"$prefix"_cpg.raw.vcf CG">> $slurms/$prefix.bissnp.sh

sbatch $slurms/"$prefix".bissnp.sh

done
