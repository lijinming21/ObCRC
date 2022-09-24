library(Boruta)
library(data.table)
library(dplyr)
library(tibble)
library(pROC)
library(caret)
library(randomForest)
library(glue)
# remotes::install_github("Tong-Chen/YSX")
library(YSX)
setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

filter.spe <- function(x) {
  return(x[grep(x, pattern = "CAG|sp\\.|unclassified|_bacterium", invert = TRUE)])
}

convert_feature <- function(feature) {
  feature[grep(feature, pattern = "X[0-9]|X\\.")] <- gsub(pattern = "X", replacement = "", feature[grep(feature, pattern = "X[0-9]|X\\.")])
  return(feature)
}

boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision),
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>%
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}

meta_crc_preparation_mupphin <- function() {
  library(MMUPHin)
  data("CRC_abd", "CRC_meta")
  
  test.mat <- read.table(paste(target[index], ".txt", sep = ""))
  test.meta <- read.table("meta.txt")
  test.meta$Cancer = ifelse(test.meta$Cancer == "Health", "CTRL", "CRC")
  
  rownames(CRC_abd) <- gsub(rownames(CRC_abd), pattern = "s__", replacement = "")
  table(rownames(test.mat) %in% rownames(CRC_abd))
  test.mat <- test.mat[filter.spe(rownames(test.mat)),]
  test.mat <- test.mat[rownames(test.mat) %in% rownames(CRC_abd),]
  test.mat <- round(test.mat/1000000, 6)
  
  dim(test.mat)
  CRC_abd <- CRC_abd[rownames(test.mat), ]
  dim(CRC_abd)
  
  CRC_meta$Cancer <- ifelse(CRC_meta$study_condition == "control", "CTRL", "CRC")
  CRC_meta$Cancer <- factor(CRC_meta$Cancer)
  CRC_meta$Group = NA
  CRC_meta <- within(CRC_meta, {
    Group[BMI < 24] <- "Normal"
    Group[BMI >= 24 & BMI < 28] <- "Overweight"
    Group[BMI >= 28] <- "Obesity"
  })
  
  table(CRC_meta$Group)
  
  table(CRC_meta$multigroup)
  CRC_meta
  
  CRC_abd <- cbind(CRC_abd, test.mat)
  
  colnames(CRC_meta)
  CRC_meta <- CRC_meta[, 28:30]
  colnames(test.meta)
  test.meta$studyID <- "JL2022"
  test.meta <- test.meta[, c(22, 1, 2)]
  CRC_meta <- rbind(CRC_meta, test.meta)
  
  fit_adjust_batch <- adjust_batch(feature_abd = CRC_abd,
                                   batch = "studyID", 
                                   covariates = "Cancer",
                                   data = CRC_meta)
  
  CRC_meta$multigroup <- NA
  CRC_meta <- within(CRC_meta, {
    multigroup[Group=="Normal" & Cancer == "CTRL"] <- "N-CTRL"
    multigroup[Group=="Overweight" & Cancer == "CTRL"] <- "Ov-CTRL"
    multigroup[Group=="Obesity" & Cancer == "CTRL"] <- "Ob-CTRL"
    multigroup[Group=="Normal" & Cancer == "CRC"] <- "N-CRC"
    multigroup[Group=="Overweight" & Cancer == "CRC"] <- "Ov-CRC"
    multigroup[Group=="Obesity" & Cancer == "CRC"] <- "Ob-CRC"
  })
  
  return(list(fit_adjust_batch$feature_abd_adj, CRC_meta))
}

crcdata <- meta_crc_preparation_mupphin()
mat <- crcdata[[1]]
meta <- crcdata[[2]]

feat.all = mat[, meta$studyID == "JL2022"]
meta.all = meta[meta$studyID == "JL2022", ]

discover <- read.table("meta.txt")
meta.all = meta.all[rownames(discover)[discover$discover == 1], ]

# abundance.threshold <- 0.0000000
# abundance = apply(feat.all, 1, mean)
# table(abundance > abundance.threshold)
# feat.all <- feat.all[abundance > abundance.threshold,]

lapply(1:3, function(i){
  # group.test <- list(c("COb", "HOb"), c("COv", "HOv"), c("CN", "HN"), c("COb", "HOb", "COv", "HOv", "CN", "HN"))
  idy <- meta.all$Group %in% group.test[i]
  table(idy)
  # change uncor
  # uncor_feat <- read.table(paste("boruta/", target[index], '/', group.test[[i]], '_uncor_feature.txt', sep=''))
  # uncor_feat <- convert_feature(uncor_feat[,1])
  # rownames(feat.all) <- gsub(pattern = '\\[|\\]|-|\\(|\\)|\\:|\\;|\\/|\\ |\\,', replacement = ".", rownames(feat.all))
  feat.sig <- t(feat.all[, rownames(meta.all)[idy]]) #uncor_feat
  
  data.sig <- cbind(feat.sig, data.frame(meta.all[rownames(feat.sig), 'Cancer']))
  colnames(data.sig)[ncol(data.sig)] <- 'Group'
  Group <- factor(data.sig$Group)
  #
  timestart<-Sys.time()
  set.seed(i + 111)
  boruta <- Boruta(x = feat.sig, y = Group, pValue = 0.001, mcAdj=T,
                   maxRuns=1000)
  timeend <- Sys.time()
  runningtime <- timeend-timestart
  print(runningtime)
  boruta
  print(table(boruta$finalDecision))
  #extract feature
  boruta.Confirmed <- data.frame(Item = getSelectedAttributes(boruta, withTentative = F), Type="Confirmed")
  feature <- boruta.Confirmed$Item
  boruta.variable.imp.use <- boruta$ImpHistory[,feature]
  feature_importance <- apply(boruta.variable.imp.use,2,mean)
  feature_importance <- data.frame(sort(feature_importance,decreasing = TRUE))
  feature <- rownames(feature_importance)
  #
  feature_importance
  write.csv(feature_importance, glue("diagnosis_model/Public_Feature/{target[index]}_{group.test[i]}_feature.csv"))
  
  ##importance
  boruta.variable.imp <- boruta.imp(boruta)
  head(boruta.variable.imp)
  #
  write.csv(boruta.variable.imp, glue("diagnosis_model/Public_Feature/{target[index]}_{group.test[i]}_feature_boruta_variable_imp.csv"))
  feature_impor_plot <- sp_boxplot(boruta.variable.imp, melted=T, xvariable = "Variable", yvariable = "Importance",
                                   legend_variable = "finalDecision", #legend_variable_order = c("shadowMax", "shadowMean", "shadowMin", "Confirmed"),
                                   legend_variable_order = c("Confirmed"), manual_color_vector = "red",
                                   xtics_angle = 90, coordinate_flip =T,
                                   outlier = TRUE, notch = FALSE)
  #
  pdf(glue("diagnosis_model/Public_Feature/{target[index]}_{group.test[i]}_boruta_feature.pdf"), useDingbats = FALSE) #,width = 6, height = 9)
  plot(feature_impor_plot)
  dev.off()
})
# })