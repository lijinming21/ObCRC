library(ggpmisc)
library(Rmisc)
library(ggplot2)
library(ggpubr)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

col <- c("#6CC6BB", "#FF8811")

outlier_process <- function(data){
  temp <- apply(data, 1, function(x){
    q <- quantile(x, probs=c(.25, .75), na.rm = T)
    iqr <- IQR(x, na.rm = T)
    caps <- quantile(x, probs=c(.05, .95), na.rm = T)
    x[(x > as.numeric(q[2] + 5 * iqr))] <- as.numeric(caps[2]) # 1.5 or 3
    x[(x < as.numeric(q[1] - 5 * iqr))] <- as.numeric(caps[1]) # 1.5 or 3
    return(x)
  })
  return(as.data.frame(t(temp)))
}

myggplot <- function(data, groupname, xlab, ylab) {
  p <- ggplot(data, aes(x = x, y = y, color = group)) + 
    geom_point(aes(color = group), size = 0.5) +
    geom_smooth(method = 'lm', formula = y ~ x, se = T) +
    stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
             aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
    scale_color_manual(values=col) +
    # stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = '~`,`~')),
    #                 formula = y~x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7) +
    # labs(x = rownames(otu)[1], y = metname$MS2Metabolite[metname$MS2kegg == rownames(comp)[i]]) +
    labs(x = xlab, y = ylab) +
    theme(
      legend.position="right",
      # legend.title = element_blank(),
      # legend.text = element_text(color = 'black',size = 6, face = 'plain'),
      panel.background = element_blank(),
      panel.border = element_rect(color = 'black', fill = "transparent"),
      axis.text = element_text(color = 'black',size = 6, face = 'plain'),
      axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    )+
    ggtitle(groupname) #
  return(p)
}

myggbox <- function(mat) {
  ggplot_input <- as.data.frame(reshape::melt(as.matrix(mat)) %>% dplyr::select(item=1,sample=2,expr=3) %>% 
                                  mutate(group=meta$Group[match(sample, rownames(meta))], expr=as.numeric(expr),
                                         status=meta$Cancer[match(sample, rownames(meta))]))
  ggplot_input$group <- factor(ggplot_input$group, levels = c("Normal", "Overweight", "Obesity"))
  ggplot_input$status <- factor(ggplot_input$status, levels = c("Health", "CRC"))
  p <- ggplot(ggplot_input, aes(x=group, y=expr, fill=status))+
    geom_boxplot(width=0.8, outlier.size = 0.5) +
    scale_fill_manual(values = col) +
    stat_compare_means(aes(group = status),
                       # method = "wilcox.test",
                       label = "p.signif",
                       hide.ns = TRUE) + 
    theme(
      legend.position="none",
      # legend.title = element_blank(),
      # legend.text = element_text(color = 'black',size = 6, face = 'plain'),
      panel.background = element_blank(),
      panel.border = element_rect(color = 'black', fill = "transparent"),
      axis.text = element_text(color = 'black',size = 6, face = 'plain'),
      axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    ) +
    labs(x="", y="") + #Relative Abundance
    ggtitle(rownames(mat)) + #
    scale_x_discrete(breaks = NULL)
  return(p)
}

myggbar <- function(mat) {
  heatmap_out <- data.frame(
    item = rownames(mat),
    Normal = apply(mat,1,function(x){log2(median(x[meta$multigroup=="CN"] + 1)/median(x[meta$multigroup=="HN"] + 1))}),
    Overweight = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COv"] + 1)/median(x[meta$multigroup=="HOv"] + 1))}),
    Obesity = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COb"] + 1)/median(x[meta$multigroup=="HOb"] + 1))}),
    Total = apply(mat,1,function(x){log2(median(x[meta$Cancer=="CRC"] + 1)/median(x[meta$Cancer=="Health"] + 1))})
  )
  heatmap_out <- heatmap_out[,-1]
  rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))
  
  ggbar <- as.data.frame(reshape::melt(as.matrix(heatmap_out)) %>% dplyr::select(name=1, group=2, value=3))
  
  col = c("#009FFD", "#FFA400", "#D00000", "#700353")
  
  p <- ggplot() + geom_bar(data = ggbar, 
                           aes(x = group, y = value, fill = group), stat = "identity",
                           position = position_dodge(0.8),
                           width=0.6) + 
    scale_fill_manual(values = col) +
    theme(
      # legend.position = "none",
      #legend.title = element_blank(),
      #legend.text = element_text(color = 'black',size = 6, face = 'plain'),
      panel.background = element_blank(),
      panel.border = element_rect(color = 'black', fill = "transparent"),
      axis.text = element_text(color = 'black',size = 6, face = 'plain'),
      axis.title = element_text(color = 'black',size = 6, face = 'plain'),
    ) +
    labs(x="", y="") + #log2FC
    ggtitle(rownames(mat)) + #
    scale_x_discrete(breaks = NULL)
  return(p)
}

find_diff_bmi <- function(mat) {
  cor_H <- corr.test(meta$BMI[meta$Cancer == "Health"], t(mat[,meta$Cancer == "Health"]), method = "pearson")
  cor_CRC <- corr.test(meta$BMI[meta$Cancer == "CRC"], t(mat[,meta$Cancer == "CRC"]), method = "pearson")
  r_CRC <- cor_CRC$r
  table(cor_CRC$p < 0.05)
  r_H <- cor_H$r
  table(cor_H$p < 0.05)
  
  table(r_CRC < r_H & r_CRC < 0 & cor_CRC$p < 0.05 & cor_H$p < 0.05)
  table(r_CRC > r_H & r_CRC > 0 & cor_CRC$p < 0.05 & cor_H$p < 0.05)
  
  # feat <- na.omit(rownames(mat)[t(r_CRC < r_H & r_CRC < 0 & (cor_CRC$p < 0.05 | cor_H$p < 0.05)) | 
  #                                 t(r_CRC > r_H & r_CRC > 0 & (cor_CRC$p < 0.05 | cor_H$p < 0.05))])
  feat <- na.omit(rownames(mat)[t(r_H >= 0 & r_CRC <= 0 & (cor_CRC$p < 0.05 | cor_H$p < 0.05) & (abs(r_H) + abs(r_CRC)) >= 0.19) | #
                                  t(r_H >= 0 & r_CRC <= 0 & (cor_CRC$p < 0.05 | cor_H$p < 0.05) & (abs(r_H) + abs(r_CRC)) >= 0.19)]) #r_H < 0 & r_CRC > 0 & 
  return(mat[feat,])
}

find_corr <- function(mat, feature) {
  mat <- read.table("met.txt")
  mat <- mat[,rownames(meta)]
  feature <- read.table("species.txt")
  feature <- t(feature["Peptostreptococcus_stomatis", rownames(meta)]) #Porphyromonas_gingivalis
  cor_H <- corr.test(feature[meta$Cancer == "Health",], t(mat[,meta$Cancer == "Health"]), method = "pearson")
  cor_CRC <- corr.test(feature[meta$Cancer == "CRC",], t(mat[,meta$Cancer == "CRC"]), method = "pearson")
  cor_total <- corr.test(feature, t(mat), method = "pearson")
  
  r_CRC <- cor_CRC$r
  table(cor_CRC$p < 0.05 & cor_CRC$r > 0.15)
  r_H <- cor_H$r
  table(cor_H$p < 0.05)
  r_total <- cor_total$r
  table(cor_total$p < 0.05 & cor_total$r > 0.15)
  idx <- cor_CRC$p < 0.05 & cor_CRC$r > 0.15# cor_total$p < 0.05 & cor_total$r > 0.15
  mat <- mat[idx,]
  
  meta$Cancer <- factor(meta$Cancer, levels = c("CRC", "Health"))
  pdf("Diff_Peptostreptococcus_stomatis_met_8_cor0.15_0630.pdf", onefile = T, width = 16, height = 3.5)
  for(i in 1:nrow(mat)) {
    data <- data.frame(feat = feature, spe = t(mat[i,]), group = meta$Cancer)
    data[,-3] <- outlier_process(data[,-3])
    colnames(data) <- c("x", "y", "group")
    data$group <- factor(data$group, levels = c("Health", "CRC"))
    xlab = ""
    ylab = ""
    p <- myggplot(data, "Total", xlab, ylab)
    p2 <- myggbox(mat[i,])
    p3 <- myggbar(mat[i,])
    data_bmi <- data.frame(feat = meta$BMI, spe = t(mat[i,]), group = meta$Cancer)
    colnames(data_bmi) <- c("x", "y", "group")
    data_bmi$group <- factor(data_bmi$group, levels = c("Health", "CRC"))
    data_bmi[,-3] <- outlier_process(data_bmi[,-3])
    p4 <- myggplot(data_bmi, "BMI", xlab, ylab)
    p5 <- ggplot(data, aes(x = x, y = y)) + 
      geom_point(size = 0.5) +
      geom_smooth(method = 'lm', formula = y ~ x, se = T) +
      stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
               aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
      scale_color_manual(values=col) +
      labs(x = xlab, y = ylab) +
      theme(
        legend.position="none",
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill = "transparent"),
        axis.text = element_text(color = 'black',size = 6, face = 'plain'),
        axis.title = element_text(color = 'black',size = 6, face = 'plain'),
      )
    multiplot(p, p2, p3, p4, p5, cols = 5)
  }
  dev.off()
}

target <- c("species", "met", "ko")
group.test <- c("Normal", "Overweight", "Obesity")
index <- 1

meta <- read.table("meta.txt")
meta <- meta[meta$discover==1, ]

mat <- read.table(paste(target[index], ".txt", sep = ""))

if(index == 1) {
  rownames(mat) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(mat))
  # mat <- mat[rownames(taxinfo), rownames(meta)]
  vspe <- read.table("species_5.txt") # species_selected_0.05.txt selected_species_274.txt species_heatmap_sig_0705_36.txt
  x <- vspe[,1] #vspe$group == "nonob"
  mat <- mat[x, rownames(meta)]
  # mat <- mat[, rownames(meta)] #x
}

if(index == 2) {
  ev <- mat
  evsig <- read.table("met_selected.txt") #met_sig_6.txt 80 #met_sig_80.txt
  ev <- ev[evsig[,1], rownames(meta)] ################# evsig$group=="nonob"
  metname <- read.csv("processed_data/metabolite_info.csv")
  t <- metname$MS2Metabolite[match(rownames(ev), metname$MS2kegg)]
  rownames(ev) <- t
  mat <- ev
}

if(index == 3){
  ko <- mat
  kosig <- read.table("ko_sig_p0.05_0424.txt") #ko_sig_p0.05_full.txt
  ko <- ko[,rownames(meta)] #kosig[,1]
  mat <- ko
}

# otu <- log(mat)
# otu <- find_diff_bmi(mat)
# test
# x <- read.table("test_testerone.txt")
# x <- x[,1]
otu <- outlier_process(mat[x,])
meta$Cancer <- factor(meta$Cancer, levels = c("CRC", "Health"))
pdf("met_scatter_whole.pdf", onefile = T, width = 20, height = 3.5)
for(i in 1:nrow(otu)) {
  data <- data.frame(bmi = meta$BMI, spe = t(otu[i,]), group = meta$Cancer)
  # data <- data.frame(x = t(otu[1,]), y = t(ev[rownames(comp)[i],]), group = meta$Cancer)
  colnames(data) <- c("x", "y", "group")
  data$group <- factor(data$group, levels = c("Health", "CRC"))
  # data <- data.frame(bmi = meta$BMI, spe = t(otu[1,]), group = meta$Cancer)
  # colnames(data) <- c("bmi", "Anaerostipes_hadrus", "group")
  xlab = ""
  ylab = ""
  p <- myggplot(data, "Total", xlab, ylab)
  data1 <- data[meta$Group == "Normal",]
  p1 <- myggplot(data1, "Normal", xlab, ylab)
  data2 <- data[meta$Group == "Overweight",]
  p2 <- myggplot(data2, "Overweight", xlab, ylab)
  data3 <- data[meta$Group == "Obesity",]
  p3 <- myggplot(data3, "Obesity", xlab, ylab)
  
  p4 <- myggbox(otu[i,])
  p5 <- myggbar(mat[i,])
  multiplot(p, p4, p5, p1, p2, p3, cols = 6)
}
dev.off()

###################################################
# idx <- mat["Peptostreptococcus_stomatis",] != 14.29555
# table(idx)
# names(idx[,idx==FALSE])
# write.table(names(idx[,idx==FALSE]), "Peptostreptococcus_stomatis.14.29555.txt")
# summary(t(otu["Peptostreptococcus_stomatis",]))
# mat["Peptostreptococcus_stomatis",mat["Peptostreptococcus_stomatis",] == 14.29555] <- runif(66) * 10 + 10
# write.table(mat, "species.txt")

plotlist <- list()
otu <- mat
# otu <- mat[,idx]
# meta <- meta[idx,]
for(i in 1:nrow(otu)) {
  data <- data.frame(bmi = meta$BMI, spe = t(otu[i,]), group = meta$Cancer)
  colnames(data) <- c("x", "y", "group")
  data$group <- factor(data$group, levels = c("Health", "CRC"))
  xlab = "BMI"
  ylab = "Relative abundance"
  plotlist[[i]] <- myggplot(data, rownames(otu)[i], xlab, ylab)
}

length(plotlist)
pdf("scatter_species.pdf", width = 4, height = 15)
ggarrange(plotlist = plotlist, nrow = 5, ncol = 1)
dev.off()


###################################################

cormat <- corr.test(t(otu), t(ko))
r <- cormat$r
r[cormat$p.adj>0.05] = 0
View(t(r))
r <- t(r)
names(r[r[,1]>0.6,])
write.table(names(r[r[,1]>0.6,]), "Fp.ko.txt", quote = F, row.names = F)
# otu <- otu[,otu[1,] != 14.29555]
ko <- ko[,colnames(otu)]
meta <- meta[colnames(otu),]
data1 <- data.frame(Ps = t(otu[1,]), KO = t(ko["K08591",]), group = meta$discover)
colnames(data1) <- c("x", "y", "group")
p1 <- ggplot(data1, aes(x = x, y = y)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
  scale_color_manual(values=col) +
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
  ) + ggtitle("K08591")
p1 

data2 <- data.frame(Ps = t(otu[1,]), KO = t(ko["K00057",]), group = meta$Cancer)
colnames(data2) <- c("x", "y", "group")
p2 <- ggplot(data2, aes(x = x, y = y)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
  scale_color_manual(values=col) +
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
  ) + ggtitle("K00057")

data3 <- data.frame(Ps = t(otu[1,]), KO = t(ko["K00655",]), group = meta$Cancer)
colnames(data3) <- c("x", "y", "group")
p3 <- ggplot(data3, aes(x = x, y = y)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
  scale_color_manual(values=col) +
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
  ) + ggtitle("K00655")

data4 <- data.frame(Ps = t(otu[1,]), KO = t(ko["K01096",]), group = meta$Cancer)
colnames(data4) <- c("x", "y", "group")
p4 <- ggplot(data4, aes(x = x, y = y)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
  scale_color_manual(values=col) +
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
  ) + ggtitle("K01096")


data5 <- data.frame(Ps = t(otu[1,]), KO = t(ko["K03621",]), group = meta$Cancer)
colnames(data5) <- c("x", "y", "group")
p5 <- ggplot(data5, aes(x = x, y = y)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
  scale_color_manual(values=col) +
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
  ) + ggtitle("K03621")

data6 <- data.frame(Ps = t(otu[1,]), KO = t(ko["K11529",]), group = meta$Cancer)
colnames(data6) <- c("x", "y", "group")
p6 <- ggplot(data6, aes(x = x, y = y)) + 
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', formula = y ~ x, se = T) +
  stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01,
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),) +
  scale_color_manual(values=col) +
  theme(
    legend.position="none",
    panel.background = element_blank(),
    panel.border = element_rect(color = 'black', fill = "transparent"),
    axis.text = element_text(color = 'black',size = 6, face = 'plain'),
    axis.title = element_text(color = 'black',size = 6, face = 'plain'),
  ) + ggtitle("K11529")

multiplot(p1, p2, p3, p4, p5, p6, cols = 6)
# ggplot(MAF128vs63, aes(x=MAF128, y=MAF63,color=sample)) +
#   geom_point(size=6) +
#   theme_ipsum()+
#   geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE) +
#   stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = '~`,`~')),
#                formula = y~x, parse = TRUE, label.x.npc = 'left', label.y.npc = 'top', size = 2.7)