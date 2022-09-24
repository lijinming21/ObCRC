rm(list = ls())
library(tidyverse)
library(dplyr)
library(glue)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(MetBrewer)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")
meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1,]
mat <- read.table("met.txt")
mat <- mat[,rownames(meta)]
metname <- read.csv("processed_data/metabolite_info.csv")
# mat$name <- rownames(mat)
mat[1:5,1:5]

s.mat <- mat

index = 2
group.test <- c("Normal", "Overweight", "Obesity")

log2fc.threshold = 0
VIP.threshold = 1
pvalue.threshold = 0.05

make_diff <- function(s.mat, i) {
  VIP <- read.table(paste("met0803/VIP_", group.test[i], ".txt", sep = ""))
  
  diff_output <- data.frame(
    feature <- rownames(s.mat),
    concentration = apply(s.mat, 1, median),
    pvalue = apply(s.mat,1,function(x){wilcox.test(unlist(x[meta$Cancer == "CRC" & meta$Group == group.test[i]]),unlist(x[meta$Cancer == "Health" & meta$Group == group.test[i]]))$p.value}),
    log2fc = apply(s.mat,1,function(x){log2(median(na.omit(x[meta$Cancer == "CRC" & meta$Group == group.test[i]]))/
                                              median(na.omit(x[meta$Cancer == "Health" & meta$Group == group.test[i]])))
      }),
    group = rep(group.test[i], nrow(s.mat))
  ) 
  
  diff_output$name <- metname$MS2Metabolite[match(rownames(s.mat), metname$MS2kegg)]
  diff_output$vip <- VIP[,2][match(rownames(diff_output), VIP[,1])]
  
  diff_output <- diff_output %>% arrange(pvalue) %>% mutate(fdr=p.adjust(pvalue,method = "BH")) %>%
    mutate(label=ifelse(log2fc > log2fc.threshold, "UP", ifelse(log2fc < (-log2fc.threshold), "DOWN", "Nosig"))) %>%
    mutate(TPplotlabel = ifelse((pvalue<pvalue.threshold & abs(log2fc)>log2fc.threshold & vip > VIP.threshold), name, NA)) %>%
    mutate(TPplotcolor = as.factor(paste(label,ifelse(is.na(TPplotlabel),0,1),sep = "")))
  
  diff_output <- diff_output[!duplicated(diff_output$TPplotlabel) | is.na(diff_output$TPplotlabel),]
  return(diff_output)
}

obdiff <- make_diff(s.mat, 3)
ovdiff <- make_diff(s.mat, 2)
ndiff <- make_diff(s.mat, 1)

# temp <- rbind(obdiff[!is.na(obdiff$TPplotlabel),], ovdiff[!is.na(ovdiff$TPplotlabel),])
df <- rbind(ndiff, ovdiff) %>% rbind(obdiff)

# write.csv(df, "met0803/vocano_fdr0.05_vip1_log2fc0.csv")

# df <- read.table("met0803/vocano_met_total_0803_pvalue_0.05_VIP_1_log2fc_0.txt")

head(df,1)
df$cluster = ifelse(df$group == "Normal", 1, ifelse(df$group == "Overweight", 2, 3))
# df$signature

# 先画背景柱，根据数据log2FC的max值,min值来确定
#根据数据中log2FC区间确定背景柱长度：

col1<-data.frame(x=c(1,2,3),
                 y=rep(max(df$log2fc), 3)) # max(t2dm$log2fc), max(crc$log2fc), max(t2dmcrc$log2fc)
col2<-data.frame(x=c(1,2,3),
                 y=rep(min(df$log2fc), 3)) # min(t2dm$log2fc), min(crc$log2fc), min(t2dmcrc$log2fc)
# 绘制背景柱
p1 <- ggplot()+
  geom_col(data = col1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = col2,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)
p1

#把散点火山图叠加到背景柱上：
p2 <- ggplot()+
  geom_col(data = col1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = col2,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = df,
              aes(x = cluster , y = log2fc, color = TPplotcolor),
              size = 1,
              width =0.4)+
  scale_color_manual(labels = c("Nosig0"="No Significance","UP1"="Significantly Up","UP0"="Up","DOWN1"="Significantly Down","DOWN0"="Down"),
                     values = c("UP0"='#FAC0AE',"DOWN0"='#9BCFF0',"UP1"='#FA2311',"DOWN1"='#6175DB',"Nosig0" ='gray'
                                # "dmCRC_Character"='#00BD9D', "dmCRC_Enhancer"='#FFC09F', "dm_CRC_risk_factor"='#52489C'
                                ))+
  labs(x="", y="log2(FoldChange)")

p2

# 添加X轴的分组色块标签：
dfcol<-data.frame(x = c(1:3),
                  y = 0,
                  label = c(1:3))
# 添加分组色块标签
dfcol$group <- group.test

# 自定义分组色块的颜色
tile_color <- met.brewer("Thomas", 3)
# 在图中镶嵌色块
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.5, #####高度
                     color = "black",
                     fill = tile_color,
                     alpha = 0.8,
                     show.legend = F)+
  geom_text(data=dfcol,
            aes(x=x, y=y, label=group),
            size =3.5,
            color ="white")
p3

log2fc.threshold <- 0.585

# library(stringr)
# idx <- by(df, df$feature....rownames.s.mat., function(x){
#             if(x$fdr[x$group == "Normal"] < 0.005 & x$fdr[x$group == "Overweight"] < 0.005)
#               return(x$feature....rownames.s.mat.[1])
#   })
# idx <- as.character(idx)
# idx <- idx[idx!="NULL"]
# idx
# 
# idy <- df[,1][df$group == "Obesity" & df$pvalue < 0.05]
# idy
# c(idx, idy)

# concen10 <- df[df$group == "Normal",]
# idz <- concen10[order(concen10$concentration, decreasing = T), 1][1:120]
idz <- concen10[concen10$concentration > 30000, 1]
df$text <- ifelse(abs(df$log2fc) > log2fc.threshold & (df[,1] %in% idz), df$TPplotlabel, NA)
# df$text <- ifelse(abs(df$log2fc) > log2fc.threshold & df$concentration, df$TPplotlabel, NA)
df$text <- ifelse(!is.na(str_extract(df$text, ".*(?=;)")), str_extract(df$text, ".*(?=;)"), df$text) # "(?<=; ).*" 

#label 调整
p4 <- p3 + geom_text_repel(
  data = df,
  aes(x = cluster, y = log2fc, label = text, color = TPplotcolor), # TextName TPplotlabel
  segment.color = 'transparent',
  max.overlaps = 50) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 1.2),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 10)
  )
p4

p4 + 

ggsave("Multiple_Volvano_p0.05_log2fc0.585_VIP1.pdf", width = 10, height = 10)
write.csv(df, "met0803/vocano_fdr0.05_vip1_log2fc0.csv")


