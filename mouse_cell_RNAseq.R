###mouse RNA-seq with IL-1B or with not
###Figure 4b
setwd("D:/HELP/LQQ/敲除数据/RNA")

devtools::install_github("junjunlab/ClusterGVis")
library("GseaVis")
library("ggplot2")
library("ClusterGVis")
library("org.Mm.eg.db")
library("clusterProfiler")
library("enrichplot")
library("DESeq2")count_GSE226983 <- read.delim('GSE226983_gene_fpkm.txt',row.names = 1)
count_GSE226983 <- count_GSE226983 %>% filter(gene_biotype == "protein_coding")
count_GSE226983=aggregate(.~gene_name,mean,data = count_GSE226983[,c(1:13)])
count_GSE226983 <- column_to_rownames(count_GSE226983,var =  "gene_name")
index <- grep("^[0-9].*",rownames(count_GSE226983),invert = T)
count_GSE226983 <- count_GSE226983[index,]
count_GSE226983 <- count_GSE226983[index,]
count_GSE226983 <- count_GSE226983[rowSums(count_GSE226983) > 0.00001,]
count_GSE226983 <- na.omit(count_GSE226983)
saveRDS(count_GSE226983,"count_GSE226983.rds")
###############diff analysis###################
count1 <- count_GSE226983[,c(1:6)]
Count_condition <- factor(c(rep("NC", 3), rep("KO", 3)),levels = c("NC",'KO'))

coldata <- data.frame(row.names=colnames(count1), Count_condition)
count1 <- as.matrix(count1)

dds <- DESeqDataSetFromMatrix(countData = round(count1),colData = coldata,design = ~Count_condition)#若batch为空则design=~label
dds <- DESeq(dds)
result <- results(dds,alpha = 0.1)
summary(result)
result <- result[order(result$padj),]
NC_VS_KO_result <- data.frame(result)
NC_VS_KO_result <- na.omit(NC_VS_KO_result)
NC_VS_KO_result$SYMBOL <- row.names(NC_VS_KO_result)
write.csv(NC_VS_KO_result,"RNA_KO_VS_NC_DESeq_result.csv")

AML_genelist_up_l = bitr(unique(row.names(NC_VS_KO_result)),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
df_all<-merge(NC_VS_KO_result,AML_genelist_up_l,by="SYMBOL",all=F)
df_all <- df_all %>% dplyr::arrange(desc(log2FoldChange))
geneList <- df_all$log2FoldChange
names(geneList) <- df_all$ENTREZID
# KEGG enrich
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               minGSSize    = 5,
               pvalueCutoff = 0.1,
               verbose      = FALSE)

# transform entrizid to gene symbol
kk2 <- setReadable(kk2,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENTREZID")
saveRDS(kk2,"RNA_KO_vs_WT_kk2.rds")
write.csv(kk2@result,"RNA_KO_vs_WT_kk2.csv")

RNA_No_il2_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/RNA/RNA_KO_vs_WT_kk2.rds")
for (i in  c(RNA_No_il2_kk2@result$ID)){
  mygene <- strsplit(RNA_No_il2_kk2@result[i,]$core_enrichment,"/")[[1]]
  gseaNb(object = RNA_No_il2_kk2,
         geneSetID = i,
         subPlot = 3,
         addGene = mygene,
         kegg = TRUE ,
         rmSegment = TRUE)
  ggsave(paste0("D:/HELP/LQQ/敲除数据/RNA/RNA_No_il2_pathways/RNA_No_IL2_kk2_",i,".pdf"))
}

################with IL-1B#########################################
count1 <- count_GSE226983[,c(7:12)]
Count_condition <- factor(c(rep("NC", 3), rep("KO", 3)),levels = c('NC',"KO"))

coldata <- data.frame(row.names=colnames(count1), Count_condition)
count1 <- as.matrix(count1)

dds <- DESeqDataSetFromMatrix(countData = round(count1),colData = coldata,design = ~Count_condition)#若batch为空则design=~label
dds <- DESeq(dds)

result <- results(dds,alpha = 0.1,contrast=c("Count_condition",'KO',"NC"))
summary(result)
result <- result[order(result$padj),]
NC_VS_KO_result <- data.frame(result)
NC_VS_KO_result <- na.omit(NC_VS_KO_result)
NC_VS_KO_result$SYMBOL <- row.names(NC_VS_KO_result)

write.csv(NC_VS_KO_result,"RNA_IL2_KO_VS_NC_DESeq_result.csv")
AML_genelist_up_l = bitr(unique(row.names(NC_VS_KO_result)),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
df_all<-merge(NC_VS_KO_result,AML_genelist_up_l,by="SYMBOL",all=F)
df_all <- df_all %>% dplyr::arrange(desc(log2FoldChange))
geneList <- df_all$log2FoldChange
names(geneList) <- df_all$ENTREZID
# KEGG enrich
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'mmu',
               minGSSize    = 5,
               pvalueCutoff = 0.1,
               verbose      = FALSE)

# transform entrizid to gene symbol
kk2 <- setReadable(kk2,
                   OrgDb = "org.Mm.eg.db",
                   keyType = "ENTREZID")
saveRDS(kk2,"RNA_IL2_kk2.rds")
write.csv(kk2@result,"RNA_IL2_KO_vs_WT_kk2.csv")

RNA_IL2_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/RNA/RNA_IL2_kk2.rds")
for (i in  c(RNA_IL2_kk2@result$ID)){
  mygene <- strsplit(RNA_IL2_kk2@result[i,]$core_enrichment,"/")[[1]]
  gseaNb(object = RNA_IL2_kk2,
         geneSetID = i,
         subPlot = 3,
         addGene = mygene,
         kegg = TRUE ,
         rmSegment = TRUE)
  ggsave(paste0("D:/HELP/LQQ/敲除数据/RNA/RNA_IL2_pathways/RNA_IL2_kk2_",i,".pdf"))
}