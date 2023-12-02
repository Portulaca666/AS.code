###mouse proteomic  with IL-1B or with not
protein_csv <- read.csv("D:/HELP/LQQ/敲除数据/protein/proteomic.csv") 
protein_csv <- protein_csv[grep("^$",protein_csv$Gene_name,invert = T),]
protein_csv <- na.omit(protein_csv)
protein_csv <- aggregate(.~Gene_name,mean,data=protein_csv)
rownames(protein_csv) <- NULL
protein_csv <- column_to_rownames(protein_csv,var = "Gene_name")
library("GseaVis")
library("ggplot2")
library("ClusterGVis")
library("org.Mm.eg.db")
library("clusterProfiler")
library("enrichplot")
library("DESeq2")
###############diff analysis###################
count1 <- protein_csv[,c(1:6)]
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
write.csv(NC_VS_KO_result,"protein_KO_VS_NC_DESeq_result.csv")

AML_genelist_up_l = bitr(unique(row.names(NC_VS_KO_result)),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
df_all<-merge(NC_VS_KO_result,AML_genelist_up_l,by="SYMBOL",all=F)#使用merge合并
# check
head(df_all)
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
write.csv(kk2@result,"protein_No_IL2_kk2.csv")
saveRDS(kk2,"protein_No_IL2_kk2.rds")

protein_No_IL2_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/protein/protein_No_IL2_kk2.rds")
show_kk2_pathways <- read.csv("D:/HELP/LQQ/敲除数据/protein/protein_No_IL2_kk2 lqq.csv")

for (i in  c(show_kk2_pathways$ID)){
  mygene <- strsplit(protein_No_IL2_kk2@result[i,]$core_enrichment,"/")[[1]]
  gseaNb(object = protein_No_IL2_kk2,
         geneSetID = i,
         subPlot = 3,
         addGene = mygene,
         kegg = TRUE ,
         rmSegment = TRUE)
  ggsave(paste0("D:/HELP/LQQ/敲除数据/protein/GSEA_PATHWAY/protein_No_IL2_kk2_",i,".pdf"))
}
################with IL-1B#########################################
count1 <- protein_csv[,c(7:12)]
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
write.csv(NC_VS_KO_result,"protein_IL2_KO_VS_NC_DESeq_result.csv")
AML_genelist_up_l = bitr(unique(row.names(NC_VS_KO_result)),fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Mm.eg.db")
df_all<-merge(NC_VS_KO_result,AML_genelist_up_l,by="SYMBOL",all=F)#使用merge合并
# check
head(df_all)
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
saveRDS(kk2,"protein_IL2_kk2.rds")
write.csv(kk2@result,"protein_IL2_kk2.csv")

protein_IL2_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/protein/protein_IL2_kk2.rds")
for (i in  c(protein_IL2_kk2@result$ID)){
  mygene <- strsplit(protein_IL2_kk2@result[i,]$core_enrichment,"/")[[1]]
  gseaNb(object = protein_IL2_kk2,
         geneSetID = i,
         subPlot = 3,
         addGene = mygene,
         kegg = TRUE ,
         rmSegment = TRUE)
  ggsave(paste0("D:/HELP/LQQ/敲除数据/protein/pathways_IL2_protein/protein_IL2_kk2_",i,".pdf"))
}
