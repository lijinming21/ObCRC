library(Maaslin2)
library(MMUPHin)
library(magrittr)
library(dplyr)
library(ggplot2)
library(vegan)
library(tidyverse)

df_input_metadata <- read.table("meta.txt")
df_input_metadata <- df_input_metadata[df_input_metadata$discover == 1, ]

target <- c("species", "met", "ko")
save_path <- c("species_maaslin2", "met_maaslin2", "ko_maaslin2")
group.test <- c("Normal", "Overweight", "Obesity")

index <- 1

df_input_data <- read.table(paste(target[index], ".txt", sep = ""))
df_input_data[1:5,1:5]

for(i in 1:3) {
  idy <- df_input_metadata$Group %in% group.test[i]
  table(idy)
  
  fit_data <- Maaslin2(
    input_data = df_input_data[,idy],
    input_metadata = df_input_metadata,
    output = paste("Result/", save_path[index], "/", group.test[i], "_output", sep = ""),
    fixed_effects = c("Cancer", "Gender", "AgeYear"), #"BMI",
    reference = c("Cancer, Health"),
    max_significance = 0.1,
    min_prevalence = 0.1,
    min_abundance = 1)
  
  maaslin2_all_results <- fit_data$results # Save results table
  maaslin2_results <- maaslin2_all_results %>% filter(metadata == 'Cancer') %>% arrange(pval) # Discard covariate associations
  maaslin2_results$qval <- p.adjust(maaslin2_results$pval, method = 'BH') # FDR correction using 'BH'
  write.table(maaslin2_results, paste("Result/", save_path[index], "/", group.test[i], "_maaslin2.txt", sep = ""),
              quote = FALSE)
}

index <- 3

df_input_data <- read.table(paste(target[index], ".txt", sep = ""))
df_input_data <- df_input_data[, rownames(df_input_metadata)]
df_input_data[1:5,1:5]
group <- c("Health", "CRC")
for(i in 1:2) {
  idy <- df_input_metadata$Cancer %in% group[i] # CRC
  table(idy)
  fit_data <- Maaslin2(
    input_data = df_input_data[,idy],
    input_metadata = df_input_metadata,
    output = paste("Result/", save_path[index], "/", group[i], "_output", sep = ""),
    fixed_effects = c("Gender", "AgeYear", "BMI"), #
    # fixed_effects = c("Group", "Gender", "AgeYear", "BMI"), #
    # reference = c("Group, Normal"),
    max_significance = 0.25,
    min_prevalence = 0.01,
    min_abundance = 1)
  
  maaslin2_all_results <- fit_data$results # Save results table
  maaslin2_results <- maaslin2_all_results %>% filter(value == 'BMI') # Discard covariate associations
  maaslin2_results$qval <- p.adjust(maaslin2_results$pval, method = 'BH') # FDR correction using 'BH'
  write.table(maaslin2_results, paste("Result/", save_path[index], "/", group[i], "_maaslin2.txt", sep = ""),
              quote = FALSE)
}


lapply(1:3, function(index){
  df_input_data <- read.table(paste(target[index], ".txt", sep = ""))
  df_input_data <- outlier_process(df_input_data)
  fit_data <- Maaslin2(
    input_data = df_input_data,
    input_metadata = df_input_metadata,
    output = paste("Result/newresult/", save_path[index], "/_sixgroup", sep = ""),
    fixed_effects = c("multigroup", "Gender", "AgeYear", "BMI"), #"BMI",
    reference = c("multigroup, HN"),
    max_significance = 0.55,
    min_prevalence = 0.1,
    min_abundance = 1)
  
  maaslin2_all_results <- fit_data$results # Save results table
  maaslin2_results <- maaslin2_all_results %>% filter(metadata == 'multigroup') %>% arrange(pval) # Discard covariate associations
  maaslin2_results$qval <- p.adjust(maaslin2_results$pval, method = 'BH') # FDR correction using 'BH'
  write.table(maaslin2_results, paste("Result/newresult/", save_path[index], "/_sixgroup", "_maaslin2.txt", sep = ""),
              quote = FALSE)
}