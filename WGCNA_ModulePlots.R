

rm(list=ls());  # empty workspace
setwd("~/A07/AML_triplets/RNA_Seq/WGCNA")

library(DESeq2)
library(dplyr)
library(reshape2)
library(ggplot2)
library(WGCNA)
library(tidyr)
library(tibble)
library(ReactomePA)
library(org.Hs.eg.db)
library(clusterProfiler)

g4reac<-readRDS("~/A07/AML_triplets/RNA_Seq/WGCNA/g4reac.rds")
res<-readRDS("~/A07/AML_triplets/RNA_Seq/DiffExp/Limma/resLimma_Rel.vs.Dia.in.Mut.rds")
a<-readRDS("~/A07/AML_triplets/RNA_Seq/DiffExp/Limma/resLimma_Rel.vs.Dia.in.Mut_named_filter.rds")

#### Reactome Pathways  ################################################################################################

#### include the DE genes (limma)
ensgEntrez<-AnnotationDbi::select(org.Hs.eg.db, keys=rownames(res),
                                  keytype="ENSEMBL", columns=c("ENSEMBL","ENTREZID") )
df<-res %>%
  rownames_to_column("ENSEMBL") %>%
  mutate(logFC=ifelse(AveExpr<1, NA, ifelse(adj.P.Val>= 0.05, NA, logFC))) %>%
  dplyr::select(ENSEMBL, logFC, adj.P.Val) %>%
  inner_join(ensgEntrez) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(lfc=mean(logFC)) # some ENTREZ genes have multiple ensembl ids

lfc <-df$lfc
names(lfc)<-df$ENTREZID

# logFCs of all genes!
df2<-res %>%
  rownames_to_column("ENSEMBL") %>%
  dplyr::select(ENSEMBL, logFC, adj.P.Val) %>%
  inner_join(ensgEntrez) %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarise(lfc=mean(logFC))

lfc2 <-df2$lfc
names(lfc2)<-df2$ENTREZID

#### skyblue ############################################################################################################
# Module genes involved in Cell cycle and DNA repair UP
grey<-g4reac %>% dplyr::filter(module=="skyblue") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()

# Reactome
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.05)
aa<-as.data.frame(x) # check for count >5
aa2<-as.data.frame(x) %>% filter(p.adjust<0.05)

par(mar=c(5,4,4,2)) #make figure margins bigger
dotplot(x, showCategory=15, font.size = 8)
emapplot(x, showCategory = 10)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2)

# GOs
go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

selected_pathways<-go_result$Description[c(1,2,5,8)]

cnetplot(go, showCategory = selected_pathways, foldChange=lfc,vertex.label.cex = 1.2)



# save the plots
# reactome
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_skyblue_reactome.png", units="cm", res= 300, width=16, height=12)
dotplot(x, showCategory=15, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_skyblue_reactome.png", units="cm", res= 300, width=20, height=20)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_skyblue2_reactome.png", units="cm", res= 300, width=30, height=30)
cnetplot(x, showCategory = 15, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()
# GO
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_skyblue_go.png", units="cm", res= 300, width=16, height=12)
dotplot(go, showCategory=15, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_skyblue_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = selected_pathways, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

# save the tables
write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_skyblue.csv")
write.csv(aa, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/pathway_result_skyblue.csv")

# list of DE genes in the modules
ListGo<-go_result %>%
  filter(Count>4)%>%
  mutate(geneID= strsplit(geneID, "/")) %>%
  unnest(geneID) %>%
  pull(geneID)%>%
  unique()

de<- a %>% filter(adj.P.Val< 0.05)

DeListGo<-intersect(ListGo, de$gene_name)
write.csv(DeListGo, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/skyblue_DEgenes_go.csv")
##########################################################################################################################

##### midnightblue ######################################################################################################
# Module Immunesystem DOWN
grey<-g4reac %>% dplyr::filter(module=="midnightblue") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.05)
aa<-as.data.frame(x) # check for count >5
aa2<-as.data.frame(x) %>% filter(p.adjust<0.05)
# only two of the pathways have count > 5

par(mar=c(5,4,4,2)) #make figure margins bigger
dotplot(x, showCategory=2, font.size = 8)
emapplot(x, showCategory = 10)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )

# GOs
go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

# filter the low counts out
go_filter <-go_result[1:15,] %>%
  filter(Count>5)
# select pathways for cnet

selected_pathways <- go_filter$Description

dotplot(go, showCategory=selected_pathways, font.size = 8)

selected_pathways2 <- selected_pathways[c(1,4,6)]

cnetplot(go, showCategory = selected_pathways2, foldChange=lfc,vertex.label.cex = 1.2)


## endocytosis pathways
#"LDLRAD3/LILRB4/LRP1/ASGR1/LILRB1/ITGB2/FCGR1A/FCGR2B/FCGR1B"
#cnetplot(go, showCategory = selected_pathways[9], foldChange=lfc,vertex.label.cex = 1.2)
# endocytosis in wrong direction. should be up in cancer


#save the plots
# reactome
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_midnightblue_reactome.png", units="cm", res= 300, width=16, height=12)
dotplot(x, showCategory=2, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_midnightblue_reactome.png", units="cm", res= 300, width=20, height=20)
cnetplot(x, showCategory = 2, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

# go
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_midnightblue_go.png", units="cm", res= 300, width=16, height=12)
dotplot(go, showCategory=selected_pathways, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_midnightblue_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = selected_pathways2, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_midnightblue.csv")
write.csv(aa, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/pathway_result_midnightblue.csv")

# list of DE genes in the modules
ListGo<-go_result %>%
  filter(Count>4)%>%
  mutate(geneID= strsplit(geneID, "/")) %>%
  unnest(geneID) %>%
  pull(geneID)%>%
  unique()

de<- a %>% filter(adj.P.Val< 0.05)

DeListGo<-intersect(ListGo, de$gene_name)
write.csv(DeListGo, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/midnightblue_DEgenes_go.csv")
######################################################################################################################

##### lightgreen module #######################################################################################
# Metabolism UP
grey<-g4reac %>% dplyr::filter(module=="lightgreen") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.1)
aa<-as.data.frame(x) # check for count >5
aa2<-as.data.frame(x) %>% filter(p.adjust<0.05) %>% # only 7 with adjP< 0.05
  mutate(GeneRatio2=as.numeric(gsub(".*/","", GeneRatio))) %>%
  mutate(GeneRatio2=Count/GeneRatio2)

p<-ggplot(aa2,
          aes(x = GeneRatio2, y = reorder(Description,Count))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  ylab(NULL)
p

ggsave("dotplot_lightgreen_reactome.pdf", plot = p, width= 200, height= 100, units= "mm", device = "pdf", 
       path = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules")

cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )
#save the plots
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_lightgreen_reactome.png", units="cm", res= 300, width=20, height=20)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

# GOs
go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

dotplot(go, showCategory=15, font.size = 8)
cnetplot(go, showCategory = 5, foldChange=lfc2,vertex.label.cex = 1.2)


selected_pathways<-go_result$Description[c(1:2)]

cnetplot(go, showCategory = selected_pathways, foldChange=lfc,vertex.label.cex = 1.2)


# macroautophagy
cnetplot(go, showCategory = go_result$Description[10], foldChange=lfc,vertex.label.cex = 1.2)
cnetplot(go, showCategory = go_result$Description[10], foldChange=lfc2,vertex.label.cex = 1.2)
# 3 de genes all downregulated.

# methylation
cnetplot(go, showCategory = go_result$Description[9], foldChange=lfc,vertex.label.cex = 1.2)
cnetplot(go, showCategory = go_result$Description[9], foldChange=lfc2,vertex.label.cex = 1.2)
# JARID2 tumor supressor
# JARID2 loss reduces H3K27me3 in gene bodies of MYCN and RUNX1T1 increasing expression

# cellular amino acid metabolic process
cnetplot(go, showCategory = go_result$Description[2], foldChange=lfc,vertex.label.cex = 1.2)
cnetplot(go, showCategory = go_result$Description[2], foldChange=lfc2,vertex.label.cex = 1.2)
# LARS1
# FARS2
# EARS2
# (WARSA is diff down, but not in modules)

# cellular amino acid metabolic process
cnetplot(go, showCategory = go_result$Description[1:2], foldChange=lfc,vertex.label.cex = 1.2)



# BCKDHB
# ALDH6A1
# MCCC1


png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_lightgreen_go.png", units="cm", res= 300, width=20, height=20)
dotplot(go, showCategory=15, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_lightgreen_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = selected_pathways, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_lightgreen.csv")
write.csv(aa2, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/pathway_result_lightgreen.csv")

# list of DE genes in the modules
ListGo<-go_result %>%
  filter(Count>4)%>%
  mutate(geneID= strsplit(geneID, "/")) %>%
  unnest(geneID) %>%
  pull(geneID)%>%
  unique()

de<- a %>% filter(adj.P.Val< 0.05)

DeListGo<-intersect(ListGo, de$gene_name)
write.csv(DeListGo, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/lightgreen_DEgenes_go.csv")

# list of De genes of the top15 modules
ListGo<-go_result[1:5,] %>%
  filter(Count>4)%>%
  mutate(geneID= strsplit(geneID, "/")) %>%
  unnest(geneID) %>%
  pull(geneID)%>%
  unique()

de<- a %>% filter(adj.P.Val< 0.05)

DeListGo<-intersect(ListGo, de$gene_name)
write.csv(DeListGo, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/lightgreen_DEgenes_goTop5.csv")
###############################################################################################################

#### lightcyan module ########################################################################################
# Nothing enriched !!
grey<-g4reac %>% dplyr::filter(module=="lightcyan") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.1)
aa<-as.data.frame(x) # check for count >5
# Nothing enriched !!

go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

dotplot(go, showCategory=15, font.size = 8)

selected_pathways <- go_result$Description[c(5,6)]
cnetplot(go, showCategory = selected_pathways, foldChange=lfc,vertex.label.cex = 1.2)

#save the plots
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_lightcyan_go.png", units="cm", res= 300, width=16, height=12)
dotplot(go, showCategory=15, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_lightcyan_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = selected_pathways, foldChange=lfc,vertex.label.cex = 1.2 )
dev.off()

write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_lightcyan.csv")
#############################################################################################################

#### Modules significant with p< 0.1 ########################################################################
#### darkturquoise module ###################################################################################
# Mainly T-cell proliferation
# Tcell and lymphocyte proliferation, leukocyte and mononuclear cells ?
# none of the GOs has DE genes!
grey<-g4reac %>% dplyr::filter(module=="darkturquoise") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.05)
aa<-as.data.frame(x) # check for count >5
# count < 5 in all pathways!

# GOs
go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

# filter the low counts out
go_filter <-go_result[1:15,] %>%
  filter(Count>5)
selected_pathways <- go_filter$Description

selected_pathways2<-go_result$Description[1:10]

dotplot(go, showCategory=selected_pathways, font.size = 8)
cnetplot(go, showCategory = selected_pathways2, foldChange=lfc,vertex.label.cex = 1.2 )
# none of the GOs has DE genes!

#save the plots
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_darkturquoise_go.png", units="cm", res= 300, width=16, height=12)
dotplot(go, showCategory=selected_pathways, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_darkturquoise_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )
dev.off()

write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_darkturquoise.csv")
write.csv(aa, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/pathway_result_darkturquoise.csv")
#############################################################################################################

#### darkgrey module ########################################################################################
# erythrocyte pathways, hemoglobin
# only 2 DE genes in reactome
# only 3 DE genes in GOs
grey<-g4reac %>% dplyr::filter(module=="darkgrey") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.05)
aa<-as.data.frame(x) # check for count >5

dotplot(x, showCategory=15, font.size = 8)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )

# GOs
go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

selected_pathways2<-go_result$Description[21:29]

dotplot(go, showCategory=14, font.size = 8)  # last one <count 5
cnetplot(go, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )
cnetplot(go, showCategory = selected_pathways2, foldChange=lfc,vertex.label.cex = 1.2 )
# RHAG, AQP1, TRIM10 only DE genes in all the GOs

#save the plots
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_darkgrey_reactome.png", units="cm", res= 300, width=16, height=12)
dotplot(x, showCategory=4, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_darkgrey_reactome.png", units="cm", res= 300, width=20, height=20)
cnetplot(x, showCategory = 4, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_darkgrey_go.png", units="cm", res= 300, width=16, height=12)
dotplot(go, showCategory=14, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_darkgrey_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_darkgrey.csv")
write.csv(aa, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/pathway_result_darkgrey.csv")
#############################################################################################################

#### brown module ##########################################################################################
# Cell Cycle
# not many DE genes
grey<-g4reac %>% dplyr::filter(module=="brown") %>% na.omit() %>% dplyr::select(ENTREZID) %>% unique()
x <- enrichPathway(gene= grey$ENTREZID, readable=T, pvalueCutoff = 0.05)
aa<-as.data.frame(x) # check for count >5

dotplot(x, showCategory=15, font.size = 8)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )

# GOs
go<- enrichGO(grey$ENTREZID,OrgDb = "org.Hs.eg.db", keyType="ENTREZID", ont = "BP", pAdjustMethod = "BH", readable = T)
go<-clusterProfiler::simplify(go) 
go_result<-go@result %>%
  filter(p.adjust<0.05)

dotplot(go, showCategory=15, font.size = 8)
cnetplot(go, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2 )


#save the plots
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_brown_reactome.png", units="cm", res= 300, width=16, height=12)
dotplot(x, showCategory=15, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_brown_reactome.png", units="cm", res= 300, width=20, height=20)
cnetplot(x, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/dotplot_brown_go.png", units="cm", res= 300, width=16, height=12)
dotplot(go, showCategory=15, font.size = 8)
dev.off()
png("~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/cnet_brown_go.png", units="cm", res= 300, width=20, height=20)
cnetplot(go, showCategory = 5, foldChange=lfc,vertex.label.cex = 1.2)
dev.off()

write.csv(go_result, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/go_result_brown.csv")
write.csv(aa, file = "~/A07/AML_triplets/RNA_Seq/WGCNA/GO_Modules/pathway_result_brown.csv")
##########################################################################################################
