library(Boruta)
library(corrplot)
library(dplyr)
library(data.table)
library(randomForest)
require(mlbench)
require(caret)
library(pROC)
library(ROCR)
library(pROC)

setwd("/Users/lijinming/Documents/datamove/colorectalsurgury/202107Metagenome_metabolome")

target <- c("species", "met", "ko")
meta <- read.table("meta.txt")
group.test <- c("Normal", "Overweight", "Obesity")

oorf <- function(mat, index, test.i, validate.j, seed) {
  set.seed(seed)
  meta.all <- meta# [meta$Group %in% group.test[test.i],]
  feat.all <- mat
  
  selected.feat <- read.csv(glue("diagnosis_model/{target[index]}_{group.test[test.i]}_feature.csv"))[,1]
  feat.sig <- t(feat.all[selected.feat, rownames(meta.all)]) #uncor_feat
  data.sig <- cbind(feat.sig, data.frame(meta.all[rownames(feat.sig), 'Cancer']))
  colnames(data.sig)[ncol(data.sig)] <- 'Group'
  # Group <- factor(data.sig$Group)
  
  data <- data.sig[meta.all$discover == 1 & meta.all$Group %in% group.test[test.i],]
  data$Group <- factor(data$Group)
  if(length(validate.j) > 1) {
    validate.data <- data.sig[meta.all$discover == 1 & meta.all$Group %in% group.test[validate.j],]
  } else {
    validate.data <- data.sig[meta.all$discover != 1 & meta.all$Group %in% group.test[validate.j],] # 
  }
  validate.data$Group <- factor(validate.data$Group)
  
  control <- trainControl(method="repeatedcv", number = 5, classProbs=T, summaryFunction=twoClassSummary)
  rf_default <- train(Group~., data=data, # as.factor(Group)
                      metric="ROC",
                      # preProcess = c("center", "scale"),
                      trControl=control)
  rf_default
  best_rf <- rf_default$bestTune
  # set.seed(index + 2021)
  fold <- createFolds(y = data$Group, k = 5, list = FALSE)
  fold
  metrics <- matrix(NA,5,12)
  colnames(metrics) <- c("Sensitivity","Specificity","Pos Pred Value","Neg Pred Value","Precision",
                         "Recall","F1","Prevalence","Detection Rate","Detection Prevalence",
                         "Balanced Accuracy","AUC")
  
  for(j in 1:5){
    fold_test <- data[fold == j,]
    fold_train <- data[fold != j,]
    print(table(fold_train$Group))
    print(table(fold_test$Group))
    Group <- fold_train$Group
    fold_fit <- randomForest(as.factor(Group)~., data=fold_train, mtry=best_rf$mtry,
                             ntree=500, importance=TRUE)
    fold_pred <- predict(fold_fit, fold_test)
    result <- confusionMatrix(factor(as.vector(fold_pred)),
                              as.factor(fold_test$Group),
                              mode = "everything", positive=levels(fold_test$Group)[1]) #######
    
    metrics[j,1:11] <- as.numeric(result$byClass)
    predicted_probs <- predict(fold_fit, fold_test, type = 'prob')
    pred <- prediction(predicted_probs[,2], fold_test$Group)
    auc <- performance(pred, 'auc')
    #print(auc@y.values[[1]])
    metrics[j,12] <- auc@y.values[[1]]
  }
  
  best_index <- which(metrics[,'AUC']==max(metrics[,'AUC']))[1]
  fold_test_best <- data[fold == best_index,]
  fold_train_best <- data[fold != best_index,]
  
  Group <- fold_train_best$Group
  best_model <- randomForest(as.factor(Group)~., data=fold_train_best,
                             mtry=best_rf$mtry,
                             ntree=500, importance=TRUE)
  best_model
  result_list <- list("model" = best_model,"metrics" = metrics,"best_fold"=best_index)
  print(result_list)
  
  test.pred.prob <- predict(best_model, fold_test_best, type='prob')
  validate.data$Group <- factor(validate.data$Group)
  
  validate.pred <- predict(best_model, validate.data)
  validate.pred.prob <- predict(best_model, validate.data, type='prob')
  
  validate.roc <- roc(validate.data$Group, validate.pred.prob[,1]) #, levels = c("H_Ob", "C_Ob"), print.auc.y = 40)
  validate.auc <- validate.roc$auc
  
  x <- list(validate.data$Group, validate.pred.prob[,1], validate.auc)
  return(x)
}

colvec <- c("#009FFD", "#FFA400", "#D00000") #c("#16CAB2", "#16697A", "#FF990A", "#FF0022")

ROC_Test_Validate <- function(x, i, start.plot = FALSE) {
  # group <- ifelse(length(testgroup)>1, 3, testgroup)
  # groupname <- ifelse(length(testgroup)>1, "TOTAL", ifelse(testgroup == 1, "CTRL", "T2DM"))
  vgroup <- x[[1]]
  vprob <- x[[2]]
  validate.roc <- roc(vgroup, vprob) #, levels = c("H_Ob", "C_Ob"), print.auc.y = 40)
  validate.auc <- validate.roc$auc
  if(start.plot) {
    plot(validate.roc, col = colvec[i], lwd = 2,
         lty = 1, xlab = "False Positive Rate", ylab = "True Positive Rate")
  } else {
    plot.roc(vgroup, vprob, legacy.axes = TRUE,
             col = colvec[i], lwd = 2, add = TRUE, lty = 1)
  }
  
  legend.line <- paste(group.test[i]," (AUC=", round(validate.auc,4),")", sep="")
  return(legend.line) # retur
}

index = 3
mat <- read.table(paste(target[index], ".txt", sep = ""))
meta <- read.table("meta.txt")
c = c(1,2,3)

lapply(1:3, function(test.i){
  n.normal <- oorf(mat, index, test.i, 1) 
  n.overweight <- oorf(mat, index, test.i, 2) 
  n.obesity <- oorf(mat, index, test.i, 3) 
  
  pdf(file = glue("diagnosis_model/ROC_Curve/{target[index]}_model_{group.test[test.i]}.pdf"))
    auc.normal <- ROC_Test_Validate(n.normal, 1, start.plot = TRUE)
    auc.overweight <- ROC_Test_Validate(n.overweight, 2)
    auc.obesity <- ROC_Test_Validate(n.obesity, 3)
    legend("bottomright", c(auc.normal, auc.overweight, auc.obesity), 
           lwd = 2, lty=1,
           col = colvec, bty="n")
  dev.off()
})

# 20 repeated Self Validation
AUC.Result = matrix(NA, nrow = 9, ncol = 20)
for(k in 1:20) {
  for(test.i in 1:3){ #  = test.meta[rownames(discover)[discover$discover == 1], ]
    seed = k
    n.normal <- oorf(mat, index, test.i, 1, seed) 
    n.overweight <- oorf(mat, index, test.i, 2, seed) 
    n.obesity <- oorf(mat, index, test.i, 3, seed) 
    AUC.Result[3 * test.i - 2, k] <- n.normal[[3]]
    AUC.Result[3 * test.i - 1, k] <- n.overweight[[3]]
    AUC.Result[3 * test.i, k] <- n.obesity[[3]]
  }
}

pdf(glue("{target[index]}_boxplot_20repeated_Self_Validation.pdf"))
auc_boxplot(AUC.Result)
dev.off()




index = 3
mat <- read.table(paste(target[index], ".txt", sep = ""))
c = c(1,2,3)
auclist = NA
pdf(file = glue("diagnosis_model/ROC_Curve/{target[index]}_model_leave_one_group_out.pdf"))
  n.normal <- oorf(mat, index, 1, c[-1]) # auto_rf(mat, index, test.i, 1)
  n.overweight <- oorf(mat, index, 2, c[-2]) 
  n.obesity <- oorf(mat, index, 3, c[-3]) 
  auclist[1] <- ROC_Test_Validate(n.normal, 1, TRUE)
  auclist[2] <- ROC_Test_Validate(n.overweight, 2)
  auclist[3] <- ROC_Test_Validate(n.obesity, 3)
  legend("bottomright", c(auclist[1], auclist[2], auclist[3]), 
         lwd = 2, lty=1,
         col = colvec, bty="n")
dev.off()


