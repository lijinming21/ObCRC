library(stringr)
library(devtools)
library(ade4)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)
library(ggsci)
library(RColorBrewer)
library(data.table)
library(tibble)
library(vegan)

alpha_index <- function(x,  tree = NULL, base = exp(1)) {
  result <- data.frame(richness = rowSums(x > 0), chao1 = estimateR(x)[3, ], ace = estimateR(x)[5, ],
                       shannon = diversity(x, index = 'shannon'), simpson = diversity(x, index = 'simpson'),
                       pielou = diversity(x, index = 'shannon') / log(estimateR(x)[1, ], base), 
                       goods_coverage = 1 - rowSums(x == 1) / rowSums(x))
  if(!is.null(tree)) {#PD_whole_tree
    pd <- pd(x, tree, include.root = FALSE)
    result <- cbind(result, pd[ ,1])
  }
  return(result)
}

mat <- read.table("species.txt")
groupinfo <- read.table("meta.txt")
groupinfo <- groupinfo[groupinfo$discover==1,]
mat <- mat[,rownames(groupinfo)]

x <- round(t(mat))
alpha.index <- data.frame(richness = rowSums(x > 0), chao1 = estimateR(x)[3, ], ace = estimateR(x)[5, ],
                          shannon = diversity(x, index = 'shannon'), simpson = diversity(x, index = 'simpson'),
                          pielou = diversity(x, index = 'shannon') / log(estimateR(x)[1, ], exp(1)), 
                          goods_coverage = 1 - rowSums(x == 1) / rowSums(x),
                          group = groupinfo$Group, cancer = groupinfo$Cancer,
                          multigroup = groupinfo$multigroup)

table(rownames(alpha.index) == colnames(mat))
manual_color_vector = c("#009FFD", "#FFA400", "#D00000")

gginput <- reshape2::melt(alpha.index, id.vars = c("group", "cancer", "multigroup"), measure.vars = c("richness", "chao1", "ace", "shannon", "simpson", "pielou", "goods_coverage"))
gginput$group <- factor(gginput$group, levels = c("Normal", "Overweight", "Obesity"))
gginput$multigroup <- factor(groupinfo$multigroup, levels = c("HN", "CN", "HOv", "COv", "HOb", "COb"))

my_comparison <- combn(as.character(unique(groupinfo$multigroup)), 2, simplify=FALSE)
unique(gginput$variable)
my_comparison <- my_comparison[c(8,12,3)]
plotlist <- list()
wiltestpvalue <- NA
comp <- NA
shannon <- gginput[gginput$variable == "shannon",]

for(i in 1:length(unique(gginput$variable))){
  i=4
  test <- gginput[gginput$variable==unique(gginput$variable)[i],]
  idj <- NA
  for(j in 1:length(my_comparison)) {
    wiltestpvalue[j] <- wilcox.test(test[test$multigroup == my_comparison[[j]][1],]$value, test[test$multigroup==my_comparison[[j]][2],]$value)$p.value
    comp[j] <- paste(my_comparison[[j]][1], my_comparison[[j]][2] , sep = "vs")
  }
  wiltestpvalue
  comp
  ptable <- data.frame(Comparison = comp, pvalue = wiltestpvalue)

  idj <- ifelse(wiltestpvalue < 0.05, 1, 0)
  
  plotlist[[i]] <- ggplot(test,aes(x=multigroup, y=value, fill=group))+
    geom_violin(trim=FALSE, aes(linetype=NA)) +
    geom_boxplot(width = 0.25, outlier.size = 0.25) +
    # geom_point(position = position_jitterdodge(),size=0.3)+
    stat_compare_means(comparisons = my_comparison[idj==1],
                       method = "t.test",
                       label = "p.signif",
                       hide.ns = TRUE)+
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.x = element_text(size = 8, angle = 45, vjust = 0.5))+
    labs(title = unique(gginput$variable)[i]) +
    scale_fill_manual(values = manual_color_vector)
}

length(plotlist)
filename <- "Alpha_index_boxplot_unadj.pdf"
pdf(filename, width = 8, height = 8)
print(ggarrange(plotlist = plotlist[1:7], nrow = 2, ncol = 4))
dev.off()

