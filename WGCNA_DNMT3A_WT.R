

# blast count needs to go in model i DE calling as well 
# that will not change the normalized counts matrix
# but when filtering out the NA genes in res it will change vsdMAt--> different modules in wgcna when design changed
# 

rm(list=ls());  # empty workspace
setwd("/home/weser/AML_triplets/RNA_Seq")

library(DESeq2)
library(dplyr)
library(reshape2)
library(ggplot2)
library(WGCNA)
library(tidyr)

anno<-readRDS("/data/htp/A07/AML_triplets/RNA_Seq/AnnoAndCounts/anno_CN.rds")
counts<-readRDS("/data/htp/A07/AML_triplets/RNA_Seq/AnnoAndCounts/counts_CN.rds")

# Patients with DNMT3a mutation only
anno<-anno[anno$DNMT3A.Status=="Wildtype",]
anno$DNMT3A.Status<-as.factor(as.character(anno$DNMT3A.Status))
counts<-counts[,colnames(counts) %in% anno$SampleName]

##### DEseq object
colnames(counts)==rownames(anno)
#anno<-anno[names(counts),]

# blast count in the model
# this only matters for the removal of NA genes
# the real DE genes will be added from limma
dds<-DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = anno, design =~ Stage + BlastCount)

dds<-estimateSizeFactors(dds) # to correct for library size
sizeFactors(dds)
dds<-estimateDispersions(dds)

# make matrix of normalized counts
countMatrix <- counts(dds, normalized=T)

design(dds)
# wald test
dds<-nbinomWaldTest(dds)

mcols(mcols(dds),use.names=T)
resultsNames(dds)

res<-DESeq2::results(dds, name = "Stage_c_vs_a")
#saveRDS(res, "/data/htp/A07/RNA_Seq/AML_triplets/DNMT3A_patients/res.rds")

# variance stabilization
vsd <- varianceStabilizingTransformation( dds )
vsdMat <- assay(vsd)

# get rid of NA genes (too small p value)
nonNAgenes<-rownames(res[is.na(res$padj)==F,])
vsdMat<-vsdMat[nonNAgenes,]

# Batch Effect Correction: Blast count
# batch effect correction
modmatrix <- model.matrix(~Stage, data=anno)
mat_noBatch<-limma::removeBatchEffect(vsdMat, covariates = anno[,c("BlastCount")], design=modmatrix)

##### WGCNA #######################################################################################################
# 1 cleaning up the data
#WGCNA wants a samples as rows!!
mat<- t(mat_noBatch)
gsg<- goodSamplesGenes(mat) # filter samples with too many missing entries, zero variance genes, weight below threshold

gsg$allOK # should be TRUE

# 2 Choose a set of soft thresholding powers (exponent in equation)
powers <- c(c(1:14), seq(from = 15, to=21, by=2))
allowWGCNAThreads(5)
## Allowing multi threading with up to 10 threads.
#sft <- pickSoftThreshold(mat,powerVector = powers )
#saveRDS(sft, "/data/htp/A07/AML_triplets/RNA_Seq/WGCNA/sft_DNMT3A_WT.rds")
# no soft threshold could be picked
# reduce the Rsquared to 0.7
sft <- pickSoftThreshold(mat,powerVector = powers, RsquaredCut = 0.7)
saveRDS(sft, "/data/htp/A07/AML_triplets/RNA_Seq/WGCNA/sft_DNMT3A_WT_0.7.rds")


# Soft Threshold Power
rsq<-ggplot(sft$fitIndices, aes(x=Power, y=(-sign(slope)*SFT.R.sq),label=Power ))+ geom_text(col="red")+
  ylab("Scale Free Topology Model \n Fit,signed Rˆ2")+
  xlab("Soft Threshold (power)")+
  ggtitle("Scale independence")+
  geom_vline(xintercept = sft$powerEstimate,col="red",alpha=0.5)+ geom_hline(yintercept = 0.85,alpha=0.5)+
  theme_minimal()
# Scale-free topology fit index as a function of the soft-thresholding power

# Mean connectivity as a function of the soft-thresholding power
df<-sft$fitIndices %>% tidyr::gather("type","Connectivity",median.k.,mean.k.,max.k.) 
conne<-ggplot(df, aes(x=Power, y=Connectivity ,label=Power, col=type ))+
  geom_point()+
  scale_y_log10()+
  xlab("Soft Threshold (power)")+ ggtitle("Connectivity")+ theme_minimal()
cowplot::plot_grid(rsq,conne,nrow = 1,labels = LETTERS[1:2],rel_widths = c(0.4,0.6))

softPower <-sft$powerEstimate
softPower

# Connectivity
k <-softConnectivity(mat, power= softPower)
names(k)<-colnames(mat)
#look at connectivity in histogram 
ggplot(data.frame(k=k),aes(x=k))+geom_histogram(fill=NA,col=1)+theme_minimal() 
#plot the fit to the scale free topology

mat<-mat[,rank(-k,ties.method="first")]

# Scale Free Plot
scaleFreePlot(k)

adjacency <- adjacency(mat, power = softPower)
saveRDS(adjacency,file="/data/htp/A07/AML_triplets/RNA_Seq/WGCNA/adjacency_DNMT3A_WT.rds") 

#calculate the dissimlaity matrix based on TOM
TOM <- TOMsimilarity(adjacency)
colnames(TOM) = colnames(mat)
rownames(TOM) = colnames(mat) 

saveRDS(TOM,file="/data/htp/A07/AML_triplets/RNA_Seq/WGCNA/tom_DNMT3A_WT.rds") 

collectGarbage() #cleaning R
