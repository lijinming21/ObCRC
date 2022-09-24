library(lattice)
library(Formula)
library(readxl)
library(ggpubr)
library(grid)
library(vcd)
library(tibble)
library(RColorBrewer)

library(caret)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(data.table)
library(devtools)

library(biomaRt)
library(ComplexHeatmap)
library(circlize)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

index <- 1

meta <- read.table("meta.txt")
meta <- meta[meta$discover == 1, ]

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("Normal", "Overweight", "Obesity")

mat <- read.table(paste(target[index], ".txt", sep = ""))
mat <- mat[, rownames(meta)]

if(index == 2){
  mat <- log2(mat+1)
}

{
ndiff <- read.table(paste("Result/", save_path[index], "/", group.test[1], "_maaslin2.txt", sep = ""), 
                    header = T, row.names = 1)
ovdiff <- read.table(paste("Result/", save_path[index], "/", group.test[2], "_maaslin2.txt", sep = ""), 
                     header = T, row.names = 1)
obdiff <- read.table(paste("Result/", save_path[index], "/", group.test[3], "_maaslin2.txt", sep = ""), 
                     header = T, row.names = 1)

num <- 6 # 6 pvalue / 8 fdr
threshold = 0.05
table(ndiff[,num] < threshold)
table(ovdiff[,num] < threshold)
table(obdiff[,num] < threshold)

ndiff <- ndiff[ndiff[,num] < threshold,]$feature
ovdiff <- ovdiff[ovdiff[,num] < threshold,]$feature
obdiff <- obdiff[obdiff[,num] < threshold,]$feature

# spe.chosen <- unique(c(ndiff, ovdiff, obdiff))

heatmap_out <- data.frame(
  item = rownames(mat),
  Normal = apply(mat,1,function(x){log2(median(x[meta$multigroup=="CN"] + 1)/median(x[meta$multigroup=="HN"] + 1))}),
  Overweight = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COv"] + 1)/median(x[meta$multigroup=="HOv"] + 1))}),
  Obesity = apply(mat,1,function(x){log2(median(x[meta$multigroup=="COb"] + 1)/median(x[meta$multigroup=="HOb"] + 1))})
)

heatmap_out <- heatmap_out[,-1]
# heatmap_out <- heatmap_out[abs(heatmap_out[,1]) != Inf & abs(heatmap_out[,2]) != Inf & abs(heatmap_out[,3]) != Inf,]
rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))

table(rownames(heatmap_out) %in% ndiff)

ndiff <- convert_feature(ndiff)
ovdiff <- convert_feature(ovdiff)
obdiff <- convert_feature(obdiff)

# feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
# rownames(heatmap_out) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(heatmap_out))
# print(feature)
# print(rownames(heatmap_out[feature,]))

heatmap_out <- na.omit(heatmap_out)

# Heatmap(heatmap_out[rownames(heatmap_out) %in% obdiff, ], name = "log2FC", 
#         show_column_names = TRUE, show_row_names = FALSE,
#         column_names_side = "top", show_column_dend = TRUE, show_row_dend = TRUE,
#         row_names_gp = gpar(fontsize = 10),
#         column_names_rot = 0, column_names_centered = TRUE,
#         # col = colorRamp2(c(-max, -max/2, 0, max/2, max), c("#109393", "#063737", "white", "#A8201A", "#A01F18"))
#         col = c("#063737", "white", "#A01F18")
#         )

sig_out <- data.frame(
  item = rownames(heatmap_out),
  Normal = ifelse(rownames(heatmap_out) %in% ndiff, TRUE, FALSE),
  Overweight = ifelse(rownames(heatmap_out) %in% ovdiff, TRUE, FALSE),
  Obesity = ifelse(rownames(heatmap_out) %in% obdiff, TRUE, FALSE)
)

rownames(sig_out) <- sig_out[,1]
sig_out <- sig_out[,-1]

# sig_heatmap <- data.frame(
#   item = rownames(heatmap_out),
#   Normal = ifelse(rownames(heatmap_out) %in% ndiff, heatmap_out$Normal, 0),
#   Overweight = ifelse(rownames(heatmap_out) %in% ovdiff, heatmap_out$Overweight, 0),
#   Obesity = ifelse(rownames(heatmap_out) %in% obdiff, heatmap_out$Obesity, 0)
# )
# 
# rownames(sig_heatmap) <- sig_heatmap$item
# sig_heatmap <- sig_heatmap[alteration$species,]
# dim(sig_heatmap)
# sig_heatmap$Signature <- alteration$group
# write.table(sig_heatmap, "species_sig_heatmap_0.05.txt", quote = F)

}
# 预处理完成，可跳下一步
##########

# library(UpSetR)
# proj <- c("species", "met", "ko")
# input <- c(
#   Normal = as.numeric(table(sig_out[,1])[2]),
#   Overweight = as.numeric(table(sig_out[,2])[2]),
#   Obesity = as.numeric(table(sig_out[,3])[2]),
#   "Normal&Overweight" = as.numeric(table(sig_out[,1] & sig_out[,2])[2]),
#   "Normal&Obesity" = as.numeric(table(sig_out[,1] & sig_out[,3])[2]),
#   "Overweight&Obesity" = as.numeric(table(sig_out[,2] & sig_out[,3])[2]),
#   "Normal&Overweight&Obesity" = as.numeric(table(sig_out[,1] & sig_out[,2] & sig_out[,3])[2])
# )
# 
# # Plot
# pdf(paste(proj[index], "_upset.pdf", sep = ""))
# upset(fromExpression(input),
#       nintersects = 40,
#       nsets = 3,
#       order.by = "freq",
#       decreasing = T,
#       mb.ratio = c(0.6, 0.4),
#       number.angles = 0,
#       text.scale = 1.1,
#       point.size = 2.8,
#       line.size = 1
# )
# dev.off()

library(ggvenn)
v1 = as.numeric(table(sig_out[,1])[2])
v2 = as.numeric(table(sig_out[,2])[2])
v3 = as.numeric(table(sig_out[,3])[2])
v12 = as.numeric(table(sig_out[,1] & sig_out[,2])[2])
v13 = as.numeric(table(sig_out[,1] & sig_out[,3])[2])
v23 = as.numeric(table(sig_out[,2] & sig_out[,3])[2])
v123 = as.numeric(table(sig_out[,1] & sig_out[,2] & sig_out[,3])[2])

H <-list('Normal'=c(v1 - v12 - v13 + v123, v12 -v123, v13 - v123, v123),
         'Overweight'=c(v2 - v12 - v23 + v123, v12 - v123, v23 - v123, v123),
         'Obesity'=c(v3 - v13 - v23 + v123, v23 - v123, v13 - v123, v123))

manual_color_vector = c("#009FFD", "#FFA400", "#D00000")

pdf(paste(target[index], "unadj_venn_p0.05.pdf", sep = ""), width = 4, height = 4)

ggvenn(H, show_elements=TRUE,
       stroke_color="black",
       fill_color = manual_color_vector,
       fill_alpha = 0.6,
       stroke_size = 1,
       set_name_size = 6,
       text_color = "black",
       text_size = 5,
       stroke_linetype="solid")

dev.off()

idx <- (heatmap_out[,1] > heatmap_out[,2] & heatmap_out[,2] > heatmap_out[,3] & abs(heatmap_out[,3])>0.5 |
          heatmap_out[,1] < heatmap_out[,2] & heatmap_out[,2] < heatmap_out[,3]  & abs(heatmap_out[,1])>0.5)
x1 <- rownames(mat)[rownames(mat) %in% rownames(sig_out[sig_out[,1] & sig_out[,2] & sig_out[,3],])] #共有 & idx
x2 <- rownames(mat)[rownames(mat) %in% rownames(sig_out[sig_out[,1] & sig_out[,2] & !sig_out[,3],])] #非肥胖 & idx
x3 <- rownames(mat)[rownames(mat) %in% rownames(sig_out[!sig_out[,1] & !sig_out[,2] & sig_out[,3],])] #仅肥胖  & idx
filter.spe <- function(x) {
  return(x[grep(x, pattern = "CAG|sp\\.|unclassified|_bacterium", invert = TRUE)])
}
x1 %>% filter.spe()
x2 %>% filter.spe()
x3 %>% filter.spe()

common <- data.frame(species = filter.spe(x1),
                     group = "common")
# common <- common[-2,]
non_ob <- data.frame(species = filter.spe(x2),
                     group = "nonob")
ob <- data.frame(species = filter.spe(x3),
                 group = "ob")
# ob <- ob[c(1:5,7,10:18),]
x4 <- rbind(common, non_ob) %>% rbind(ob)

write.table(x4, "species_selected_p0.05.txt")

# heatmap
s.mat <- apply(mat, 1, scale)
s.mat <- t(s.mat)
colnames(s.mat) <- colnames(mat)

# table(rownames(meta) == colnames(mat))
# meta$Location <- ifelse(meta$Location == "/", NA, meta$Location)
# col_2 <- colorRamp2(c(0, 1), c("#2FBF71", "#ED7D3A"))
column_ha = HeatmapAnnotation(Status = meta$Cancer,
                              BMI = factor(meta$Group, levels = c("Normal", "Overweight", "Obesity")),
                              Age = factor(meta$Age, levels = c("Young", "Elder")),
                              Gender = factor(meta$Gender, levels = c("Male", "Female")),
                              # Location = factor(meta$Location, levels = c("Left", "Right", "Rectal")),
                              # Differentiation = factor(meta$Differentiation, levels = c("Well-Moderate", "Poor")),
                              # Lymphatic_Metastasis = meta$Lymphatic_Invasion,
                              # Vascular_Invasion = meta$Vascular_Invasion,
                              # Stage = factor(meta$T),
                              # MSI = meta$MS,
                              col = list(Status = c("CRC" = "#1F1F1F", "Health" = "#F0F3BD"),
                                         BMI = c("Normal" = "#009FFD", "Overweight" = "#FFA400", "Obesity" = "#D00000"),
                                         Age = c("Young" = "#1F1F1F", "Elder" = "#F0F3BD"),
                                         Gender = c("Male" = "#1F1F1F", "Female" = "#F0F3BD")
                                         # Location = c("Left" = "#C64191", "Right" = "#62BFED", "Rectal" = "#297373"),
                                         # Differentiation = c("Well-Moderate" = "#1F1F1F", "Poor" = "#F0F3BD"),
                                         # Lymphatic_Metastasis = c("Positive" = "#1F1F1F", "Negative" = "#F0F3BD"),
                                         # Vascular_Invasion = c("Positive" = "#1F1F1F", "Negative" = "#F0F3BD"),
                                         # Stage = c("0" = "#F0F3BD", "1" = "#DEBAC0", "2" = "#C64191", "3" = "#62BFED",
                                         #           "4" = "#297373")
                                         # MSI = c("MSI-H" = "#1F1F1F", "MSS" = "#F0F3BD")
                              )
)

# col = list(Status = c("CRC" = "#028090", "Health" = "#F0F3BD"),
#            BMI = c("Normal" = "#009FFD", "Overweight" = "#FFA400", "Obesity" = "#D00000"),
#            Age = c("Young" = "#028090", "Elder" = "#F0F3BD"),
#            Gender = c("Male" = "#028090", "Female" = "#F0F3BD"),
#            # Location = c("Left" = "#C64191", "Right" = "#62BFED", "Rectal" = "#297373"),
#            Differentiation = c("Well-Moderate" = "#028090", "Poor" = "#F0F3BD"),
#            Lymphatic_Metastasis = c("Positive" = "#028090", "Negative" = "#F0F3BD"),
#            # Vascular_Invasion = c("Positive" = "#028090", "Negative" = "#F0F3BD"),
#            Stage = c("0" = "#F0F3BD", "1" = "#DEBAC0", "2" = "#C64191", "3" = "#62BFED",
#                      "4" = "#297373"),
#            MSI = c("MSI-H" = "#028090", "MSS" = "#F0F3BD")

signature <- read.table("met_selected.txt") # species_selected_top.txt
row_ha = rowAnnotation(bar = signature$group ,
                       col = list(bar = c("common" = "#B3CBB9", "nonob" = "#84A9C0", "ob" = "#6A66A3"))
)

sig <- convert_feature(signature[,1])
p <- Heatmap(s.mat[sig,], name = "Z-score", border = TRUE, #rownames(s.mat) %in% 
        show_column_names = FALSE, show_row_names = TRUE,
        top_annotation = column_ha,
        left_annotation = row_ha,
        # column_split = factor(meta$multigroup, levels = c("HN", "HOv", "HOb", "CN", "COv", "COb"), 
        #                       labels = c("N-CTRL", "Ov-CTRL", "Ob-CTRL", "N-CRC", "Ov-CRC", "Ob-CRC")),
        column_split = factor(meta$multigroup, levels = c("HN", "CN", "HOv", "COv", "HOb", "COb"), 
                              labels = c("N-CTRL", "N-CRC", "Ov-CTRL", "Ov-CRC", "Ob-CTRL", "Ob-CRC")),
        cluster_column_slices = FALSE,
        column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 10),
        row_order = sig,
        column_names_rot = 0, column_names_centered = TRUE,
        #col = c("#15099A", "white", "#B7151B")
        col = colorRamp2(c(-2,0,2), c("#0F8B8D", "white", "#A8201A"))
)
p

pdf(paste(target[index], "_whole_heatmap.pdf", sep = ""), width = 20, height = 10) # ,
print(p)
dev.off()

## diff_heatmap all significant and a increase/decrease trend or opposite fold change
table(sig_out[,1] & sig_out[,2] & sig_out[,3]) #sig_out[,1] & sig_out[,2] & sig_out[,3]
idx1 <- sig_out[,1] & sig_out[,2] & sig_out[,3]
(idx2 <- which(idx1 & (heatmap_out[,1] < heatmap_out[,2]) & (heatmap_out[,2] < heatmap_out[,3])))
(idx3 <- which(idx1 & (heatmap_out[,1] > heatmap_out[,2]) & (heatmap_out[,2] > heatmap_out[,3])))
(idx4 <- which(idx1 & (heatmap_out[,1] * heatmap_out[,3] < 0))) # & (abs(heatmap_out[,1]) + abs(heatmap_out[,3]) > 0.5)
idx <- unique(c(idx2, idx3, idx4))
rownames(sig_out[unique(idx),])

# significant and log2foldchange > 1 in either group
table( (sig_out[,3] & abs(heatmap_out[,3]) > 1) )
table( (sig_out[,2] & abs(heatmap_out[,2]) > 1) )
table( (sig_out[,1] & abs(heatmap_out[,1]) > 1) )

table((sig_out[,1] & abs(heatmap_out[,1]) > 1) | (sig_out[,2] & abs(heatmap_out[,2]) > 1) | (sig_out[,3] & abs(heatmap_out[,3]) > 1))
(idy <- which((sig_out[,1] & abs(heatmap_out[,1]) > 1) | (sig_out[,2] & abs(heatmap_out[,2]) > 1) | (sig_out[,3] & abs(heatmap_out[,3]) > 1)))

rownames(sig_out[unique(idy),])

(idy1 <- which(sig_out[,3] & abs(heatmap_out[,3]) > 1))
x <- rownames(sig_out[unique(idy1),])
#

x <- rownames(sig_out[unique(c(idy, idx)),])

table(sig_out[,1] | sig_out[,2] | sig_out[,3])

idx <- (sig_out[,1] & sig_out[,2] & sig_out[,3])
# heatmap <- heatmap_out[idx,]
# sig_out <- sig_out[idx,]
# idx <- which(sig_out[,3])

# (idx1 <- which(!sig_out[,1] & !sig_out[,2] & sig_out[,3])) # 
# (idsig <- sig_out[,3])
# (idx2 <- which(idsig & (heatmap_out[,1] * heatmap_out[,3] < 0))) #  & abs(heatmap_out[,3]) > 1)
# (idx3 <- which(idsig & heatmap_out[,1] < heatmap_out[,2] & heatmap_out[,2] < heatmap_out[,3]))
# (idx4 <- which(idsig & heatmap_out[,1] > heatmap_out[,2] & heatmap_out[,2] > heatmap_out[,3])) #  & abs(heatmap_out[,3]) > 1
# idx <- unique(c(idx1, idx2, idx3, idx4))

# different trend
(idx1 <- which(sig_out[,3] & ((heatmap_out[,1] * heatmap_out[,3] < 0)))) #  & abs(heatmap_out[,3]) > 0.2
(idsig <- sig_out[,1] | sig_out[,2] | sig_out[,3])
(idx2 <- which(idsig & (heatmap_out[,1] * heatmap_out[,3] < 0))) # & (abs(heatmap_out[,1]) + abs(heatmap_out[,3]) > 0.5)
(idx3 <- which(sig_out[,3] & (heatmap_out[,1] < heatmap_out[,2]) & (heatmap_out[,2] < heatmap_out[,3]))) #  & abs(heatmap_out[,3]) > 0.2
(idx4 <- which(sig_out[,3] & (heatmap_out[,1]) > (heatmap_out[,2]) & (heatmap_out[,2]) > (heatmap_out[,3]))) #  & abs(heatmap_out[,3]) > 0.2
idx <- unique(c(idx1, idx2, idx3, idx4))
# idx <- unique(c(idx1, idx2))

unique(idx)
signif <- sig_out[unique(idx),] #1
heatmap <- heatmap_out[unique(idx),] #1

write.table(rownames(heatmap) ,paste(target[index], "_heatmap_sig_0705_36.txt", sep = ""))
pdf(paste(target[index], "_heatmap_sig_0705.pdf", sep = ""), height = 8, width = 4)
print(
  Heatmap(heatmap, name = "log2FC",
          # left_annotation = ha,
          show_column_names = TRUE,# show_row_names = FALSE,
          column_names_side = "top", show_column_dend = FALSE, show_row_dend = FALSE,
          row_names_gp = gpar(fontsize = 7),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 0, column_names_centered = TRUE,
          # col = colorRamp2(c(-max, -7, 0, 7, max), c("#063737", "#109393", "white", "#A8201A", "#A01F18")),
          col = colorRamp2(c(-2,0,2), c("#0F8B8D", "white", "#A8201A")),
          column_order  = c("Normal", "Overweight", "Obesity"),
          row_order = rownames(heatmap), #rownames(heatmap_out[unique(idx),]),
          # width = unit(10, "cm"), height = unit(40, "cm"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(signif[i,j]) {
              sig <- ifelse(heatmap[i,j] > 0, "+", "-")
              grid.text(sig, x, y, gp = gpar(fontsize = 5))
            }
          }
  )
)
dev.off()
