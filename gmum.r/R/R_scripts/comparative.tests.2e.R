#I HATE R

# General results: SVMLight lack of control over convergence time, by default much longer convergence in RBF case

# Tested: cv time, final model time, accuracy, #sv, cv scores similarity

library(e1071)
library(doMC)
registerDoMC() # parallel foreach
library(caret)
library(SparseM)
library(mlbench)
library(knitr)
#source("caret.svm.models.R")


# Download libsvm datasets to data_local folder:
# australian     diabetes   german.numer  sonar_scale
# breast-cancer  fourclass  heart         splice
# Read binary datasets
datasets = lapply(list.files(path="data_local"), function(x){ ds=read.matrix.csr(paste("data_local",x,sep="/")); ds$name=x; ds})

get.results <- function(name, data, class.column.name, params, model, seed=777, ...){
  
  class.column.index <- which(colnames(data)==class.column.name)
  p = 0.66
  inTraining <- 1:round(nrow(data) * p)
  training <- data[ inTraining,]
  testing  <- data[-inTraining,]
  
  params <- params$C
  models <- sapply(params, function(x) SVM(formula(paste(class.column.name,"~ .", sep=" ")), prep="2e", training, kernel="linear", C=x, ...))
  
  l <-c()
  row_length <- 0
  for(model in models){
    class.column.index < which(colnames(data)==class.column.name)
    pred_test <- predict(model, subset(testing, select=-c(class.column.index)))
    pred_train <- predict(model, subset(training, select=-c(class.column.index)))
    cf_test <- confusionMatrix(pred_test, testing[, class.column.index])
    cf_train <- confusionMatrix(pred_train, training[, class.column.index])
    
    score_row = c(name, as.numeric(cf_train$overall[1]), as.numeric(cf_test$overall[1]), model$SV, model$C, model$iterations)
    row_length <- length(score_row)
    l <- c(l, score_row)
  }
  result <- data.frame(matrix(l, nrow=length(params), ncol=row_length, byrow=TRUE))
  colnames(result) <- c("ExperimentName", "trAcc", "tstAcc", "SVs", "C", "iterations")
  result
}

# Coarse CV grid
get.all.results <- function(ds, fit.svmlight=TRUE, fit.small.grid=TRUE){
  subset <- c(1)
  
  C <- 10^(-6:7)
  
  if(fit.small.grid){
    C <- 10^(-1:1)
  }
  
  model.names=list("gmum.r::svm.linear")
  model.calls=list(caret.gmumSvmLinear)
  model.tuneGrids = list(expand.grid(C=C))
  
  model.names=model.names[subset]
  model.calls=model.calls[subset]
  model.tuneGrids = model.tuneGrids[subset]
  
  df <- as.data.frame(as.matrix(ds$x))
  df$y = ds$y
  df$name
  
  res <- matrix(, nrow = 0, ncol = 6)
  for(model.index in 1:length(model.names)){
    name <- model.names[[model.index]]
    if(grepl("^gmum", name)){
      name.experiment <- paste("libsvm_",ds$name, sep="")
      res <- rbind(res, get.results(name.experiment, df, "y", model.tuneGrids[[model.index]], model.calls[[model.index]], verbosity=0, lib='libsvm'))
      if(fit.svmlight){
        name.experiment <- paste("svmlight_",ds$name, sep="")
        res <- rbind(res, get.results(name.experiment, df, "y", model.tuneGrids[[model.index]], model.calls[[model.index]], verbosity=0, lib='svmlight'))
      }
    }
  }
  
  print(kable(res, format = "markdown"))
  
}

for(ds in datasets){
  get.all.results(ds, fit.svmlight = TRUE, fit.small.grid = TRUE)
}
