library(ade4)
library(vegan)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(dplyr)
require(readxl)
library(ropls)

rundate <- "plsda0803_nolog"
mat <- read.table("met.txt")
mat[1:5,1:5]
dim(mat)

groupinfo <- read.table("meta.txt")
mat <- mat[,groupinfo$discover == 1]
groupinfo <- groupinfo[groupinfo$discover == 1,]
s.mat <- as.data.frame(apply(mat, 1, function(x){scale(x)}))
rownames(s.mat) <- colnames(mat)
s.mat[1:5,1:5]

my_comparison <- combn(as.character(unique(groupinfo$multigroup)), 2, simplify=FALSE)
my_comparison <- my_comparison[-c(4,5,7,9,10,11)]
my_comparison

mat <- as.data.frame(s.mat)

for(i in 1:length(my_comparison){
  if(i %in% c(3,5,6)) {
    my_comparison[[i]]
    idx <- groupinfo$multigroup %in% my_comparison[[i]]
    met.plsda <- opls(s.mat[idx,], groupinfo[idx,]$multigroup, predI = 1, orthoI = 1, permI = 100,
                      crossvalI = 7)
    # met.plsda@summaryDF
    met.plsda@modelDF
    #O-PLSDA
    sample.score <- met.plsda@scoreMN %>% cbind(met.plsda@orthoScoreMN[,1]) %>%
                      as.data.frame() %>%
                      mutate(group = groupinfo[idx,]$multigroup)#, o1 = met.plsda@orthoScoreMN[,1])
    # OPLSDA
    # sample.score <- met.plsda@orthoScoreMN %>%
    #   as.data.frame() %>%
    #   mutate(group = groupinfo[idx,]$multigroup)#, o1 = met.plsda@orthoScoreMN[,1])
    
    colnames(sample.score) <- c("p1", "p2", "group")
    ggplot(sample.score, aes(p1, p2, color = group)) +
      geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) +
      geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
      geom_point() +
      labs(title = paste0("p=",met.plsda@summaryDF$pQ2),
           x = paste('t1(', round(met.plsda@modelDF$R2X[1]*100), '%)', sep = ''),
           y = paste('o1(', round(met.plsda@modelDF$R2X[2]*100), '%)', sep = '')) + 
      stat_ellipse(level = 0.95, linetype = 'solid', 
                   size = 1, show.legend = FALSE) +
      scale_color_manual(values = c('#84DCC6','#EB8884')) +
      theme_bw() +
      theme(
            legend.position = "none", # right
            legend.title = element_blank(),
            legend.text = element_text(color = 'black',size = 15, face = 'plain'),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_text(color = 'black',size = 15, face = 'plain'),
            axis.title = element_text(color = 'black',size = 15, face = 'plain'),
            axis.ticks = element_line(color = 'black')
      )
    
    ggsave(paste(rundate, my_comparison[[i]][1], "vs", my_comparison[[i]][2],".pdf", sep=""), width=5,height=5)
    out_met <- met.plsda@vipVn[order(met.plsda@vipVn, decreasing = T)]
    
    write.table(sample.score, paste(rundate, my_comparison[[i]][1], "vs", my_comparison[[i]][2],"_OPLSDA.txt", sep=""))
    write.table(out_met, paste(rundate, my_comparison[[i]][1], "vs", my_comparison[[i]][2],"VIP.txt", sep=""),
                quote = FALSE, col.names = FALSE)
  }
}
