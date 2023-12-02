##figure 1d
library(tibble)
library(dplyr)
library(pheatmap)
setwd("D:/HELP/LQQ/LQQ/GEO")
normal_normalized <- read.delim("D:/HELP/LQQ/LQQ/GEO/GSE114007_normal_normalized.counts.txt")
OA_normalized <- read.delim("D:/HELP/LQQ/LQQ/GEO/GSE114007_OA_normalized.counts.txt")
normalized_data <- dplyr::inner_join(normal_normalized,OA_normalized,by = "symbol") %>% dplyr::select(-c("Average.Normal","Max.x","Average.OA","Max.y")) %>% column_to_rownames(.,var = "symbol")
normalized_data <- normalized_data[,grep("Cart",colnames(normalized_data),invert = T)]
data <- normalized_data[c("DDX5","MMP13","ADAMTS5","SOX9"),]
data$gene=rownames(data)
data=melt(data,id="gene")
colnames(data) <- c("gene","sample","value")
data$group <- gsub("[0-9.*$]","",data$sample)
data$group <- gsub("\\_","",data$group)
library(ggdist)
library(ggpubr)
violin_gene <- read.csv("D:/HELP/LQQ/LQQ/火山图基因.csv",header = F)
violin_gene <- violin_gene$V1
#########差异分析###############
sample_info <- read.csv("D:/HELP/LQQ/LQQ/sample.csv",header = T)
sample_info$combat <- c(rep("cart",18),rep("normal",20))
SraRunTable <- read.csv("D:/HELP/LQQ/LQQ/SraRunTable.txt")
sample_info$Instrument <- SraRunTable$Instrument[match(sample_info$sampleID,SraRunTable$Run)]
sample_info$name <- sample$sample[match(sample_info$SRR,sample$sampleID)]
sample_info <- na.omit(sample_info)
rownames(sample_info) <- sample_info$name
sample_info$combat <- c(rep("Normal",10),rep("OA",10))

##bulkRNA_data was processed by fastq --> bam --> count --> TPM
library(limma)
library(edgeR)
DGElist <- DGEList(counts = bulkRNA_data, group = factor(sample_info$combat))
keep_gene <- rowSums(cpm(DGElist) > 1) >= 2
DGElist <- DGElist[keep_gene,keep.lib.sizes = FALSE]
DGElist <- calcNormFactors(DGElist)
design <- model.matrix( ~0 + factor(sample_info$combat))
rownames(design) <- colnames(DGElist)
colnames(design) <- levels(factor(sample_info$combat))
v <- voom(DGElist, design, plot = TRUE, normalize = "quantile")
fit <- lmFit(v, design)
cont.matrix <- makeContrasts(contrasts = c('OA-Normal'), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
nrDEG_limma_voom = topTable(fit2, coef = 'OA-Normal', n = Inf)
nrDEG_limma_voom = na.omit(nrDEG_limma_voom)
head(nrDEG_limma_voom)
##plot
library(ggplot2)
library(ggrepel)
nrDEG <- nrDEG_limma_voom
nrDEG$change <- ifelse(nrDEG$adj.P.Val < 0.05 & abs(nrDEG$logFC) > 1,
                       ifelse(nrDEG$logFC > 1,'UP','DOWN'),
                       'NOT')
violin_gene <-  violin_gene[violin_gene %in% rownames(nrDEG)]
nrDEG$sign <- rownames(nrDEG)
nrDEG <- nrDEG %>% filter(sign %in% violin_gene)
nrDEG$sign <- rownames(nrDEG)
nrDEG <- na.omit(nrDEG)
table(nrDEG$sign)
library(ggrepel)
ggplot(data= nrDEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(alpha=0.4,aes(size = abs(logFC))) +
  theme_bw(base_size = 15) +
  theme(plot.title=element_text(hjust=0.5), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + 
  geom_hline(yintercept=1.30103,linetype=4) +
  geom_vline(xintercept=c(-1,1) ,linetype=4 ) +
  scale_color_manual(name = "",
                     values = c('brown','steelblue','gray'),
                     limits = c("UP", "DOWN", "NOT")) +
  geom_label_repel(aes(label=sign),
                   fontface="bold",
                   color="grey50",
                   box.padding=unit(0.35, "lines"),  
                   point.padding=unit(0.5, "lines"), 
                   segment.colour = "grey50", 
                   force = T) +
  labs(title = 'OA DEG volcano')
ggsave("OA_DEG_volcano.pdf",width = 15,height = 10)

##figure 1e
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
compaired <- list(c("normal","OA"))
g <- ggplot(data, aes(x = group, y = value, color = group, fill = group)) +
  scale_y_continuous(breaks = 1:9) +
  scale_color_manual(values = my_pal, guide = "none") +
  scale_fill_manual(values = my_pal, guide = "none") + theme_bw()
g + 
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)
  ) + geom_signif(comparisons = compaired,map_signif_level = F,test = "wilcox.test") +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 3
  )  +  facet_wrap(~gene,nrow = 2,scales = "free_y")  -> p1
print(p1)
ggsave("4genes_expr.pdf")