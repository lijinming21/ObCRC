# BiocManager::install("preprocessCore")
# 加载包
library(WGCNA)
library(igraph)
library(psych)
library(impute)
library(corrplot)

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

group.test <- c("Normal", "Overweight", "Obesity")

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")

metname <- read.csv("processed_data/metabolite_info.csv")

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]

otu <- read.table("species.txt")
ev <- read.table("met.txt")
ko <- read.table("ko.txt")

otusig <- read.table("species_10.txt") #species_selected_p0.05.txt species_heatmap_sig_0705_36.txt selected_species.txt 48 # selected_species_274.txt # spe_sig.txt # species_heatmap_sig_0518.txt # species_heatmap_sig_0524_simple.txt
evsig <- read.table("met_selected.txt") # met_sig_18_0708.txt met_sig_p0.05_0514.txt 80 met_sig_6.txt
kosig <- read.table("ko_sig_anno_full.txt") #ko_sig_p0.05_full.txt ko_sig_p0.05_0424.txt ko_sig_anno_full_0706.txt
# kosig <- read.table("ko_database.txt", header = T, row.names = 1)
kosig <- kosig[kosig$ko %in% rownames(ko),]

# rownames(kosig) <- 1:nrow(kosig)
# x <- table(kosig$type)
# x <- x[order(x, decreasing = T)]
# x
# table(kosig$type %in% names(x[c(1,2,4,7,8,9,13)]))
# kosig <- kosig[kosig$type %in% names(x[c(1,2,4,7,8,9,13)]),]
# kosig <- kosig[!duplicated(kosig$ko),]
# dim(kosig)

table(rownames(otu) %in% otusig[,1])
rownames(otu) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(otu))
table(rownames(otu) %in% otusig[,1])

colnames(otusig)[1] <- "feature"
colnames(evsig)[1] <- "feature"
colnames(kosig)[1] <- "feature"
otusig$group2 <- "species" # paste0("species_", otusig$group)
evsig$group2 <- "metabolites" # paste0("metabolites_", evsig$group)
kosig$group <- "ko"
kosig$group2 <- kosig$subtype # type
kosig <- subset(kosig, select = c("feature", "group", "group2"))
# kosig <- kosig[c(1:113, 119:124),]

otu <- read.table("species_10_mat.txt")
otu <- otu[otusig[1:8,1], rownames(meta)]
ev <- ev[evsig[,1], rownames(meta)]
ko <- ko[kosig[,1], rownames(meta)]
# otu <- otu[na.omit(match(taxinfo$species, rownames(otu))),]
t <- metname$MS2Metabolite[match(rownames(ev), metname$MS2kegg)]
rownames(ev) <- t

evsig$feature <- metname$MS2Metabolite[match(evsig[,1], metname$MS2kegg)]
feature.temp <- rbind(otusig, evsig) %>% rbind(kosig)

temp <- rbind(otu, ev)
temp <- rbind(temp, ko)

# start
rvalue <- 0.2
ko.rvalue <- 0.6
pvalue <- 0.0005
point.del <- 1
withname <- TRUE
pdf(paste("correlation/", "r", rvalue, "_rko", ko.rvalue, "_p", pvalue, "_del", point.del, "_", withname,
          "_0909.pdf", sep = ""),
    width = 10, height = 28)
par(mfrow=c(3, 1)) #, mar=c(1,1,1,1)

layouts = c("layout_with_fr", "layout_with_kk", "layout_with_dh",
            "layout_with_gem", "layout_as_star", "layout_as_tree",
            "layout_in_circle", "layout_on_grid")

lapply(1:3, function(i){
  OTU <- t(temp)
  idy <- meta$Group %in% group.test[i]
  table(idy)
  OTU <- OTU[idy,]

  occor = corAndPvalue(OTU)
  
  occor.r = occor$cor # 取相关性矩阵R值
  occor.p = occor$p # 取相关性矩阵p值
  
  table(occor.p < pvalue & abs(occor.r) > rvalue)
  occor.r[occor.p > pvalue | abs(occor.r) < rvalue] = 0 
  
  idx1 <- 1:nrow(otu)
  idx2 <- (nrow(otu)+1):(nrow(otu)+nrow(ev)) #ev
  idx3 <- (nrow(otu)+nrow(ev)+1):nrow(temp) #ko
  idx4 <- (nrow(otu)+1):nrow(temp) #ev+ko
  # occor.r[idx1, idx1] = 0
  occor.r[idx2, idx2] = 0
  occor.r[idx3, idx3] = 0
  occor.r[idx4, idx4] = 0
  #################################################
  for(k in idx1) {
    for(j in idx3) {
      occor.r[k, j] = ifelse(abs(occor.r[k, j]) < ko.rvalue, 0, occor.r[k, j])
      occor.r[j, k] = ifelse(abs(occor.r[j, k]) < ko.rvalue, 0, occor.r[j, k])
    }
  }

  igraph = graph_from_adjacency_matrix(occor.r, mode = "undirected", weighted=TRUE, diag=FALSE)

  if (point.del > 0) {
    bad.vs = V(igraph)[igraph::degree(igraph) < point.del]
    igraph = igraph::delete.vertices(igraph, bad.vs)
  }
  igraph.weight = E(igraph)$weight
  E(igraph)$weight = NA
  set.seed(123)
  E.color = igraph.weight
  E.color = ifelse(E.color>0, "#D04539",ifelse(E.color<0, "#3B94BA", "grey"))
  E(igraph)$color = as.character(E.color)
  V(igraph)$size = 5
  
  # set vertices color
  rownames(feature.temp) <- feature.temp$feature
  igraph.col = feature.temp[V(igraph)$name,]
  
  colist <- c("#EBF2FA", "#00A878", "#D8F1A0", "#F3C178", "#FE5E41", "#1E91D6", "#086788", "#5B5941", "#E9D2F4", "#340068", "#D88C9A")
  if(length(colist) <= length(unique(igraph.col$group2))) {
    colist <- rep(colist, length.out = length(unique(igraph.col$group2)))
  } else {
    colist <- colist[1:length(unique(igraph.col$group2))]
  }
  names(colist) <- unique(igraph.col$group2)
  ver.col <- factor(igraph.col$group2, levels = names(colist), labels = colist)
  V(igraph)$color = as.character(ver.col)
  set.seed(123)
  
  l=do.call(layouts[4], list(igraph))
  if(withname) {
    plot(igraph, main=paste(group.test[i], " Co-occurrence network", sep = ""),
       vertex.frame.color="black", edge.arrow.size = 1, # vertex.label = NA,
       edge.lty=1, edge.curved=FALSE, vertex.label.dist = 0)
  } else {
    plot(igraph, main=paste(group.test[i], " Co-occurrence network", sep = ""),
         vertex.frame.color="black", vertex.label = NA, edge.arrow.size = 1,
         edge.lty=1, edge.curved=FALSE)
  }
  legend("right", legend = names(colist), fill = colist)
  # x <- as_data_frame(igraph, what="edges")
  # y <- as_data_frame(igraph, what="vertices")
})

dev.off()

