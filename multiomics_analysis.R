###multiomics  proteomic and RNA-seq analysis
#############NO IL2####################
RNA_KO_WT_DESeq <- read.csv("D:/HELP/LQQ/敲除数据/RNA/RNA_KO_VS_NC_DESeq_result.csv")
protein_KO_WT_DESeq <- read.csv("D:/HELP/LQQ/敲除数据/protein/protein_KO_VS_NC_DESeq_result.csv")
link_KO_WT_DESeq <- RNA_KO_WT_DESeq %>%  inner_join(protein_KO_WT_DESeq, by = "SYMBOL")
##x =  RNA ； y = protein
write.csv(link_KO_WT_DESeq,"link_KO_WT_DESeq.csv")
library(ggpubr)
gggrepl_gene <- read.csv("D:\\HELP\\LQQ\\敲除数据\\双组学联合\\link_KO_WT_DESeqlqq.csv")
gggrepl_gene$X <- NULL
ggplot(link_KO_WT_DESeq)+ geom_point(aes(log2FoldChange.x,log2FoldChange.y),shape = 21, color="grey", size=5,fill = "#FFF0F5")+
  geom_point(data = gggrepl_gene,aes(log2FoldChange.x,log2FoldChange.y),shape = 21, color="grey", size=5,fill = "red")+
  geom_label_repel(data = gggrepl_gene,aes(log2FoldChange.x, log2FoldChange.y, fill="blue", 
                                           label=X.x), fontface="bold", color="white", 
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+ theme_classic(base_size = 16)+
  geom_hline(yintercept=0,linetype=2)+    
  geom_vline(xintercept=0,linetype=2)+
  xlab("RNA-seq log2FoldChange")+
  ylab("Protein-seq log2FoldChange")+
  ggtitle("RNA-seq and protein Diff analysis")+
  NoLegend()
ggsave("No_IL2_link.pdf")
############# With   IL2####################
RNA_KO_WT_DESeq <- read.csv("D:/HELP/LQQ/敲除数据/RNA/RNA_IL2_KO_VS_NC_DESeq_result.csv")
protein_KO_WT_DESeq <- read.csv("D:/HELP/LQQ/敲除数据/protein/protein_IL2_KO_VS_NC_DESeq_result.csv")
link_KO_WT_DESeq <- RNA_KO_WT_DESeq %>%  inner_join(protein_KO_WT_DESeq, by = "SYMBOL")
##x =  RNA ； y = protein
write.csv(link_KO_WT_DESeq,"IL2_link_KO_WT_DESeq.csv")
library(ggpubr)
gggrepl_gene <- read.csv("D:\\HELP\\LQQ\\敲除数据\\双组学联合\\IL2_link_KO_WT_DESeq lqq.csv")
gggrepl_gene$X <- NULL

ggplot(link_KO_WT_DESeq)+ geom_point(aes(log2FoldChange.x,log2FoldChange.y),shape = 21, color="grey", size=5,fill = "#FFF0F5")+
  geom_point(data = gggrepl_gene,aes(log2FoldChange.x,log2FoldChange.y),shape = 21, color="grey", size=5,fill = "red")+
  geom_label_repel(data = gggrepl_gene,aes(log2FoldChange.x, log2FoldChange.y, fill="blue", 
                                           label=X.x), fontface="bold", color="white", 
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+ theme_classic(base_size = 16)+
  geom_hline(yintercept=0,linetype=2)+  
  geom_vline(xintercept=0,linetype=2)+
  xlab("RNA-seq log2FoldChange")+
  ylab("Protein-seq log2FoldChange")+
  ggtitle("RNA-seq and protein Diff analysis")+
  NoLegend()

ggsave("IL2_link.pdf")
###############pathways  link#######################################
protein_IL2_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/protein/protein_IL2_kk2.rds")
RNA_IL2_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/RNA/RNA_IL2_kk2.rds")

link_il2_DESeq <- protein_IL2_kk2@result %>%  inner_join(RNA_IL2_kk2@result, by = "Description")
ggplot(link_il2_DESeq)+ geom_point(aes(enrichmentScore.x,enrichmentScore.y),shape = 21, color="grey", size=5,fill = "#FFF0F5")+
  geom_point(data = link_il2_DESeq,aes(enrichmentScore.x,enrichmentScore.y),shape = 21, color="grey", size=5,fill = "red")+
  geom_label_repel(data = link_il2_DESeq,aes(enrichmentScore.x, enrichmentScore.y, fill="blue", 
                                             label=Description), fontface="bold", color="white", 
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+ theme_classic(base_size = 16)+
  geom_hline(yintercept=0,linetype=2)+     
  geom_vline(xintercept=0,linetype=2)+
  scale_x_continuous(limits = c(0.3,1)) +
  scale_y_continuous(limits = c(0.3,1)) +
  xlab("protein_IL2 enrichmentScore")+
  ylab("RNA-IL2 enrichmentScore")+
  ggtitle("IL2-RNA-seq and protein pathway enrichment")+
  NoLegend()
ggsave("Pathways_IL2_link.pdf")


protein_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/protein/protein_No_IL2_kk2.rds")
RNA_kk2 <- readRDS("D:/HELP/LQQ/敲除数据/RNA/RNA_KO_vs_WT_kk2.rds")

link_DESeq <- protein_kk2@result %>%  inner_join(RNA_kk2@result, by = "Description")
ggplot(link_DESeq)+ geom_point(aes(enrichmentScore.x,enrichmentScore.y),shape = 21, color="grey", size=5,fill = "#FFF0F5")+
  geom_point(data = link_DESeq,aes(enrichmentScore.x,enrichmentScore.y),shape = 21, color="grey", size=5,fill = "red")+
  geom_label_repel(data = link_DESeq,aes(enrichmentScore.x, enrichmentScore.y, fill="blue", 
                                         label=Description), fontface="bold", color="white", 
                   box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50")+ theme_classic(base_size = 16)+
  geom_hline(yintercept=0,linetype=2)+  
  geom_vline(xintercept=0,linetype=2)+
  # scale_x_continuous(limits = c(0.3,1)) +
  scale_y_continuous(limits = c(0.3,1)) +
  xlab("protein enrichmentScore")+
  ylab("RNA enrichmentScore")+
  ggtitle("RNA-seq and protein pathway enrichment")+
  NoLegend()
ggsave("Pathways_No_IL2_link.pdf")
