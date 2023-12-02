setwd("D:/HELP/LQQ/rMATs_new/human/output")
files <- list.files()
for (i in c(1:35)[-c(5:29)]){
  MATs <- read.table(paste0("./",files[i]),row.names = 1,header = T)
  #MATs <- MATs %>% dplyr::filter(PValue < 0.05)
  write.table(MATs,paste0("D:/HELP/LQQ/rMATs_new/human/filtered/",files[i]),sep = "\t",quote=F, col.names = T)
}
####过滤
df <- data.frame(row.names =  c("A5SS.MATS.JC","A3SS.MATS.JC","SE.MATS.JC","RI.MATS.JC","MXE.MATS.JC"))
for (i in c("A5SS.MATS.JC","A3SS.MATS.JC","SE.MATS.JC","RI.MATS.JC","MXE.MATS.JC")){
  A3SS.MATS.JC <- read.table(paste0("./",i,".txt"),sep = "\t",header = T)
  A3SS.MATS.JC <-  A3SS.MATS.JC %>% dplyr::filter(FDR < 0.05)
  Normal_events_number <- length(A3SS.MATS.JC[A3SS.MATS.JC$IncLevelDifference > 0,]$GeneID)
  OA_events_number <- length(A3SS.MATS.JC[A3SS.MATS.JC$IncLevelDifference < 0,]$GeneID)
  df[i,1] <- Normal_events_number
  df[i,2] <- OA_events_number
}

colnames(df) <-  c("Normal","OA")
df$state <- rownames(df)
df <- reshape2::melt(df, id="state", variable.name="events", value.name = "number_events")
mycol= brewer.pal(n = 12, name = "Set3")
df$number_events <- ifelse(df$events == "OA",-df$number_events,df$number_events)
df$hjust <- ifelse(df$events == "OA",0,1)
df$state <- as.character(gsub(".MATS.JC","",df$state))
df <- df %>% dplyr::arrange(number_events)
write.csv(df,"df.csv")

##figure 1b
p1 <- ggplot(df,aes(y = state)) +
  geom_col(aes(x = number_events,fill = events),show.legend = F) +
  scale_fill_manual(values = c("#1B65B1","#E87586")) +
  scale_y_discrete(position = "right",
                   name = "types of events") +
  xlab("events numbers of OA and Normal") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.position = "left") + geom_text(aes(x= number_events,label = number_events))

pdf("D:/HELP/LQQ/rMATs_new/human/filtered/filtered_SEs.pdf",8,4)
print(p1)
dev.off()

##figure 1c
need  <- colnames(rMATs)
rMATs <- data.frame()
rMATS_list <- list()

for (i in c("A5SS.MATS.JC","A3SS.MATS.JC","SE.MATS.JC","RI.MATS.JC","MXE.MATS.JC")){
  A3SS.MATS.JC <- read.delim(paste0("./",i,".txt"))
  A3SS.MATS.JC <-  A3SS.MATS.JC %>% dplyr::filter(FDR < 0.05)
  A3SS.MATS.JC <-  A3SS.MATS.JC %>% dplyr::filter(abs(IncLevelDifference) > 0) 
  A3SS.MATS.JC$state <- gsub(".MATS.JC","",i)
  A3SS.MATS.JC <-  A3SS.MATS.JC %>% dplyr::select(GeneID,FDR,IncLevel1,IncLevel2,IncLevelDifference,state) 
  rMATs <- rbind(rMATs,A3SS.MATS.JC)
  rMATS_list[[i]] <- A3SS.MATS.JC
}


rMATs$sample <- ifelse(rMATs$IncLevelDifference>0,"Normal","OA")
ECM_genes <- read.csv("D:/HELP/LQQ/rMATs_new/human/ECM.csv",header = F)
ECM_genes <- c(ECM_genes$V1,ECM_genes$V2,ECM_genes$V3,ECM_genes$V4,ECM_genes$V5)[1:164]

rMATs <- rMATs %>% dplyr::filter(GeneID %in% ECM_genes)
rownames(rMATs) <- paste(rMATs$GeneID,rMATs$state,sep = "_",rownames(rMATs))
write.csv(rMATs,"D:/HELP/LQQ/rMATs_new/human/filtered/rMATs_out.csv")

rMATs <- read.csv("D:/HELP/LQQ/rMATs_new/human/filtered/rMATs_out.csv",row.names = 1)
rMATs_graph_data <- rMATs %>% dplyr::rename(from=GeneID,to=state,group = sample,IncLevelDifference=IncLevelDifference)
library(ggraph)
library(igraph)
rMATs_graph_data <- rMATs_graph_data %>% dplyr::arrange(desc(group),IncLevelDifference)
rMATs_graph_data <- dplyr::filter(rMATs_graph_data,from %in% ECM_genes)
rMATs_graph_data$to <- as.character(rMATs_graph_data$to)
rMATs_graph_data$group <- factor(rMATs_graph_data$group)
rMATs_graph_data <- rMATs_graph_data %>% dplyr::select(from,to,group,FDR,IncLevelDifference)
actors <- data.frame(name = c(unique(rMATs_graph_data$from),unique(rMATs_graph_data$to)),attribute = c(rep("gene",length(unique(rMATs_graph_data$from))),rep("state",4))) 
actors$label <- actors$name
#intersect_genes <- c("PLOD2","FN1")
#actors$label[c(match(intersect_genes,toupper(actors$name)),(nrow(actors)-4):nrow(actors))] <- actors$name[c(match(intersect_genes,toupper(actors$name)),(nrow(actors)-4):nrow(actors))]
graph <- graph_from_data_frame(rMATs_graph_data,directed=TRUE,vertices = actors)
p2 <- ggraph(graph, layout = 'kk') + 
  geom_edge_link(mapping = aes(edge_color = group,edge_width = abs(IncLevelDifference),), 
                 alpha = 0.3,
                 start_cap = circle(3, 'mm'), 
                 end_cap = circle(2, 'mm')) + 
  geom_node_point(size = 5,
                  mapping = aes(colour = factor(attribute)),
                  alpha = 0.5) + theme_graph()+
  scale_edge_width_continuous(range = c(0.5,2)) +
  scale_edge_color_manual(values = c("#0b6bcb","#f6798b"))+
  scale_color_manual(values = c("#f6798b","#FFA500"))+
  geom_node_text(mapping = aes(label=label),repel = T,size = 4)

library(showtext)
font_add('Arial','/Library/Fonts/Arial.ttf') 
showtext_auto() 

pdf("D:/HELP/LQQ/rMATs_new/human/filtered/human_ggraph_intersect.pdf",10,8)
print(p2)
dev.off()
saveRDS(rMATs_graph_data,"D:/HELP/LQQ/rMATs_new/human/filtered/human_rMATs_graph_data.rds")

