library(glue)
library(psych)
library(reshape2)
library(WGCNA)
library(ggrepel)
library(ggClusterNet)
library(sna)
library(network)
library(dplyr)
library(scales)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("COb", "HOb", "COv", "HOv", "CN", "HN")

metname <- read.csv("processed_data/metabolite_info.csv")

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]
otu <- read.table("species.txt")

otusig <- read.table("species_10.txt") # species_heatmap_sig_0705_36.txt selected_species.txt 48 # selected_species_274.txt # spe_sig.txt # species_heatmap_sig_0518.txt # species_heatmap_sig_0524_simple.txt

table(rownames(otu) %in% otusig[,1])
rownames(otu) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(otu))
table(rownames(otu) %in% otusig[,1])

otu <- otu[otusig[,1], rownames(meta)]

# start
rvalue <- 0.3
pvalue <- 0.05
plotlist <- list()
group.test <- c("Obesity", "Overweight", "Normal")

for(i in 1:3) {
  # idy <- meta$multigroup %in% group.test[i]
  idy <- meta$Group %in% group.test[i]
  table(idy)
  mat <- t(otu[,idy])
  mat <- na.omit(mat)
  
  occor <- corr.test(mat, mat,
                     use="pairwise",
                     method="spearman", # 可选pearson/kendall
                     adjust="fdr",
                     alpha=0.05)
  
  r_matrix <- occor$r
  p_matrix <- occor$p.adj
  
  write.csv(r_matrix, glue("{group.test[i]}_species_co_occurrence_r_matrix.csv"))
  write.csv(p_matrix, glue("{group.test[i]}_species_co_occurrence_p_adj_matrix.csv"))
  
  table(p_matrix < pvalue & abs(r_matrix) > rvalue)
  idx <- p_matrix < pvalue & abs(r_matrix) > rvalue
  r_matrix[!idx]=0
  
  netClu = data.frame(ID = colnames(mat),
                      #group = otusig$group
                      group = "species"
                      )
  netClu$group = as.factor(netClu$group)
  
  set.seed(12)
  
  result2 = PolygonClusterG(cor = r_matrix, nodeGroup = netClu, zoom = 0.8, zoom2 = 0.8) #PolygonRrClusterG
  node = result2[[1]]
  head(node)
  # ---node节点注释
  # nodes = nodeadd(plotcord = node, otu_table = otu_table, tax_table = tax_table)
  nodes <- node
  nodes$group <- otusig$group # [match(node$elements, netClu$ID)]
  # nodes$type <- feature.otu$group2[match(node$elements, feature.otu$feature)]
  if(length(group.test) == 3) {
    if(i == 1){
      foldchange = apply(otu,1,function(x){log2(median(na.omit(x[meta$multigroup == "COb"]))/(median(na.omit(x[meta$multigroup == "HOb"]))))})
    } else if (i == 2) {
      foldchange = apply(otu,1,function(x){log2(median(na.omit(x[meta$multigroup == "COv"]))/(median(na.omit(x[meta$multigroup == "HOv"]))))})
    } else if (i == 3) {
      foldchange = apply(otu,1,function(x){log2(median(na.omit(x[meta$multigroup == "CN"]))/(median(na.omit(x[meta$multigroup == "HN"]))))})
    }
  } else if(length(group.test) == 6) {
    if(i <= 2){
      foldchange = apply(otu,1,function(x){log2(median(na.omit(x[meta$multigroup == "COb"]))/(median(na.omit(x[meta$multigroup == "HOb"]))))})
    } else if (i <= 4) {
      foldchange = apply(otu,1,function(x){log2(median(na.omit(x[meta$multigroup == "COv"]))/(median(na.omit(x[meta$multigroup == "HOv"]))))})
    } else if (i <= 6) {
      foldchange = apply(otu,1,function(x){log2(median(na.omit(x[meta$multigroup == "CN"]))/(median(na.omit(x[meta$multigroup == "HN"]))))})
    }
  }
  # names(foldchange[rownames(nodes)]) == rownames(nodes)
  # foldchange <- foldchange[rownames(nodes)]
  nodes$Alteration <- ifelse(foldchange > 0, "Enriched in CRC", "Depleted in CRC")
  # nodes$updown <- ifelse(foldchange > 0, 21, 23)
  abundance = log10(apply(otu,1,median))
  nodes$abund <- abundance
  
  # colnames(nodes)
  #-----计算边
  edge = edgeBuild(cor = r_matrix, node = node)
  head(edge)
  # edge$weight <- abseedge
  # colnames(edge)[8] = "cor"
  write.csv(edge, glue("{group.test[i]}_species_co_occurrence.csv"))
  colvec <- c("#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6")
  colist <- rep(colvec, length.out = length(unique(nodes$group)))
  names(colist) <- unique(nodes$group)
  
  plotlist[[i]] <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = weight,
                                               size = 1), data = edge, alpha = 0.3) + #size = (weight/3)
    #scale_colour_manual(values = c("+" = "#D00000", "-" = "#009FFD")) +
    scale_colour_gradient2(
      low = "#009FFD",
      mid = "white",
      high = "#D00000",
      limits=c(-0.8, 0.8)
    ) +
    geom_point(aes(X1, X2, fill = Alteration, size = abund), pch = 21, data = nodes) + #, size = mean , 
    scale_fill_manual(values = c("Enriched in CRC" = "#D00000", "Depleted in CRC" = "#009FFD")) +
    guides(fill=guide_legend(override.aes=list(shape=21))) +
    # geom_text_repel(aes(X1, X2, label = elements), size = 2, nudge_y = -0.2, data = nodes,
    #                 max.overlaps = 20) +
    geom_text(aes(X1, X2, label = elements), size = 4, nudge_y = -0.08, data = nodes) +
    # scale_fill_manual(values = c("Increased" = "#D00000", "Decreased" = "#009FFD")) +
    # scale_size(range = c(4, 14)) +
    scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme(legend.position = "none") +
    theme(panel.background = element_blank()) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    theme(legend.background = element_rect(colour = NA)) +
    theme(panel.background = element_rect(fill = "white",  colour = NA)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    ggtitle(group.test[i])
}
library(ggpubr)
length(plotlist)
pdf("species_co_nolegend.pdf", width = 14, height = 5)
print(ggarrange(plotlist = plotlist, nrow = 1, ncol = 3))
dev.off()

