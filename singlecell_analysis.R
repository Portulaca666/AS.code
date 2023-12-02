######singlecell RNA-seq analysis#######
library(Seurat)
library(ggplot2)
library(harmony)
library(stringr)
lf <- list.files()
sceList <- lapply(lf,FUN = function(x){
  TPM <- read.delim(x,header = TRUE,row.names = 1)
  seu <- CreateSeuratObject(counts = log(TPM+1))
  return(seu)})

sceBig <- sceList[[1]]
for (i in 2:34){
  sceBig <- merge(sceBig,sceList[[i]])
}

sce <- SCTransform(sceBig)
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))
AS.obj = sce %>% RunHarmony("orig.ident", 
                            plot_convergence = TRUE)

ElbowPlot(sce,ndims=25)

AS.obj <- RunUMAP(AS.obj,reduction = "harmony",dims = 1:4) 
AS.obj <- FindNeighbors(AS.obj,reduction = "harmony",dims = 1:4) 
AS.obj<- FindClusters(AS.obj,resolution = 0.49)

stage <- str_split(rownames(AS.obj@meta.data),pattern = "\\_",simplify = T)[,2]
AS.obj$stage <- str_split(stage,pattern = "\\.",simplify = T)[,1]
DBseu$anno <- NA
DBseu$anno[DBseu$seurat_clusters == "0"] <- "preHTC"
DBseu$anno[DBseu$seurat_clusters == "1"] <- "EC"
DBseu$anno[DBseu$seurat_clusters == "2"] <- "HTC"
DBseu$anno[DBseu$seurat_clusters == "3"] <- "ZC"
DBseu$anno[DBseu$seurat_clusters == "4"] <- "RegC"
DBseu$anno[DBseu$seurat_clusters == "5"] <- "HomC"
DBseu$anno[DBseu$seurat_clusters == "6"] <- "FC"
DBseu$anno[DBseu$seurat_clusters == "7"] <- "ProC"
DimPlot(DBseu,group.by = "anno")
allcolour=c("#549D87",
            "#603990",
            "#86CEE9",
            "#C02A1D",
            "#898EB4",
            "#F46D43",
            "#A2D2C3",
            "#F5AE9F")
##Extended Data Fig. 1.a
DimPlot(AS.obj,group.by = "anno")+ 
  scale_fill_manual(values = allcolour) +
  scale_color_manual(values = allcolour) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(size = 10),  
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=10), 
        legend.key.size=unit(0.52,'cm') ) + 
  guides(color = guide_legend(override.aes = list(size=3.5)))

features <- c("ADAMTS5","TGFBI","TNFAIP6", #preHTC
              "C2orf82","CLEC3A","CYTL1", #EC
              "SPP1","COL10A1","BHLHE41", #HTC
              "ZNF585B","SLC5A12","HSD3BP4", #ZC 
              "CHI3L2","CRTAC1","PTGS2","MMP3","CHI3L1", #RegC
              "JUN","FOS","FOSB", # HomC
              "COL1A1","COL14A1","MMP13","COL1A2","S100A4","COL3A1","COL5A1", #FC
              "COL2A1","FGF1","NGF","KRT17","KRT16","P3H2","ERCC1"#ProC
) 

df <- data.frame(gene = features,
                 type = c(rep("preHTC",3),
                          rep("EC",3),
                          rep("HTC",3),
                          rep("ZC",3),
                          rep("RegC",5),
                          rep("HomC",3),
                          rep("FC",7),
                          rep("ProC",7)))

gene <- merge(df,p$data,by.x = "gene",by.y = "features.plot")

gene$type <- factor(gene$type ,levels = c("preHTC","EC","HTC","ZC","RegC",
                                          "HomC","FC","ProC"))
##Extended Data Fig. 1.b
ggplot(gene,aes(x=gene,y =  as.numeric(id),size = pct.exp, color = avg.exp.scaled))+
  geom_point(alpha=0.45) + 
  scale_size("pct.exp", range = c(0,6)) + 
  scale_color_gradientn(colours = viridis::viridis(50,direction = -1),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "avg.exp.scaled") +
  cowplot::theme_cowplot() + 
  ylab("group") + 
  theme_bw() +
  scale_y_continuous(breaks = 1:length(levels(gene$id)),labels = levels(gene$id),sec.axis = dup_axis())+
  facet_grid(~type, scales="free_x",space = "free")+theme_classic() +
  theme(
    axis.text.x = element_text(size=10, angle=90, hjust=1, color="black"), 
    axis.text.y = element_text(size=10, color="black"),
    axis.ticks.y = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y = element_text(size=10,colour = 'black',vjust = 2.5,hjust = 0.5), 
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2), 
    panel.spacing=unit(0.5, "mm"),
    strip.text.x = element_text(size=10,color = "black",
                                vjust = 0.5,margin = margin(b = 3,t = 3)),
    strip.background = element_rect(colour="grey30", fill="white",size = 0.3)
  )
##Extended Data Fig. 1.c
AS.obj$seurat_name = AS.obj$anno
Idents(AS.obj) <- AS.obj$seurat_name
AS.obj.seurat_name.markers <- FindAllMarkers(AS.obj, only.pos = FALSE,
                                             min.pct = 0.25,assay = 'SCT',slot = "data",
                                             logfc.threshold = 0)
# top 5 genes
top5pos <- AS.obj.seurat_name.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top5negtive <- AS.obj.seurat_name.markers %>% group_by(cluster) %>% top_n(n = -10, wt = avg_log2FC)
top10 <- rbind(top5pos,top5negtive)
DDX5_gene <- AS.obj.seurat_name.markers %>% filter(gene %in% "COL2A1")
top5pos <- rbind(top5pos,DDX5_gene)
# plot
pdf("AS.obj.seurat_name.markers.pdf",13,7)
ggplot(AS.obj.seurat_name.markers,
       aes(x = pct.1 - pct.2,y = avg_log2FC)) +
  geom_point(color = 'grey80') +
  geom_hline(yintercept = c(-0.25,0.25),lty = 'dashed',size = 1,color = 'grey50') +
  geom_text_repel(data = top5pos,
                  aes(x = pct.2 - pct.1,y = avg_log2FC,
                      label = gene,color = cluster),
                  show.legend = F,direction = 'y',
                  hjust = 1,
                  nudge_y = 0.25, 
                  force = 5, 
                  nudge_x = 0.8 - (top5pos$pct.2 - top5pos$pct.1)) +
  geom_text_repel(data = top5negtive,
                  aes(x = pct.1 - pct.2,y = avg_log2FC,
                      label = gene,color = cluster),
                  show.legend = F,direction = 'y',
                  hjust = 0,
                  force = 2.5, 
                  nudge_x = -0.8 - (top5negtive$pct.2 - top5negtive$pct.1)) +
  geom_point(data = top5pos,show.legend = F,
             aes(x = pct.1 - pct.2,y = avg_log2FC,color = cluster)) +
  scale_color_npg(name = '') +
  # x y breaks label
  scale_y_continuous(limits = c(-6,10),breaks = seq(-6,10,2)) +
  scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5)) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        strip.background = element_rect(color = NA,fill = 'grey90')) +
  xlab(expression(Delta~'Percentage Diffrence')) +
  ylab('Log2-Fold Change') +
  facet_wrap(~cluster,nrow = 1,scales = 'fixed')
dev.off()

##Extended Data Fig. 2.
butu_genes <- read.csv("D:/HELP/LQQ/LQQ/LQQ_AS_scRNAseq/补图1/butu1.csv",col.names = F)
butu_genes <- butu_genes$FALSE.
TWO_obj$seurat_name <- factor(TWO_obj$seurat_name,levels = c("ProC", "FC"))
TWO_obj$DDX5 <- TWO_obj@assays[["SCT"]]@data["DDX5",]
splots <- list()
for (i in c(toupper(butu_genes)[-19],"MMP13","ADAMTS5","ADAMTS4","SOX9","DDX17")){
  TWO_obj@meta.data$COL1A1 <- TWO_obj@assays[["SCT"]]@data[i,]
  data<-TWO_obj@meta.data
  my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
  g <- ggplot(data, aes(x = seurat_name, y = COL1A1, color = seurat_name, fill = seurat_name)) +
    scale_y_continuous(breaks = 1:9) +
    scale_color_manual(values = my_pal, guide = "none") +
    scale_fill_manual(values = my_pal, guide = "none") 
  splots[[i]] <- g + 
    geom_boxplot(
      width = .2, fill = "white",
      size = 1.5, outlier.shape = NA
    ) +
    ggdist::stat_halfeye(
      adjust = .33, ## bandwidth
      width = .67, 
      color = NA, ## remove slab interval
      position = position_nudge(x = .15)
    ) + theme_bw() + 
    geom_signif(comparisons = list(c("ProC", "FC")),y_position=c(7)) + 
    gghalves::geom_half_point(
      side = "l", 
      range_scale = .3, 
      alpha = .5, size = 3
    ) + ylab(i) + theme(panel.grid=element_blank())
}

library(patchwork)
gg <- splots[[1]]+splots[[2]]+splots[[3]]+splots[[4]]+splots[[5]]+ splots[[6]]+splots[[7]]+ plot_layout(widths = c(1, 1))
pdf("gene.pdf",5,5)
print(splots)
dev.off()
##figure 1f
genes <- c("COL1A1","COL1A2","COL2A1","DDX5")
count <- TWO_obj@assays$SCT@data[genes,]
count <- data.frame(t(count))
library(corrplot)
pdf("corrplot_square.pdf",4,4)
corrplot(cor(count), method = 'square',order = 'AOE', addCoef.col = 'black', tl.pos = 'd',
         cl.pos = 'r',col =  rev(COL2('PiYG', 10)))
dev.off()