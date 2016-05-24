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
source("caret.svm.models.R")


# Download libsvm datasets to data_local folder:
# australian     diabetes   german.numer  sonar_scale
# breast-cancer  fourclass  heart         splice
# Read binary datasets
datasets = lapply(list.files(path="data_local"), function(x){ ds=read.matrix.csr(paste("data_local",x,sep="/")); ds$name=x; ds})

get.results <- function(name, data, class.column.name, params, model, seed=777, ...){
  set.seed(seed)
  
  class.column.index <- which(colnames(data)==class.column.name)
  inTraining <- createDataPartition(data[, class.column.index], p = .75, list = FALSE)
  training <- data[ inTraining,]
  testing  <- data[-inTraining,]
  
  
  
  fitControl <- trainControl(method = "cv",
                             ## 10-fold CV...
                             number = 7,
                             ## repeated 2 times
                             repeats = 1,
                             verboseIter=FALSE
  )
  
  
  
  model <- train(formula(paste(class.column.name,"~ .", sep=" ")), data = training,
                 method = model,
                 preProc = c("center", "scale"),
                 tuneLength = 8,
                 tuneGrid = params,               
                 trControl = fitControl,
                 ...)
  
  
  class.column.index < which(colnames(data)==class.column.name)
  pred <- predict(model$finalModel, predict(model$preProcess, subset(testing, select=-c(class.column.index))))
  cf <- confusionMatrix(pred, testing[, class.column.index])
  
  list(name=name,
       trainAcc=max(model$results$Accuracy), 
       trainAccStd=model$results$AccuracySD[which.max(model$results$Accuracy)], 
       testAcc=as.numeric(cf$overall[1]),
       trainTime=as.numeric(model$times$everything[1]),
       finalModelTime=as.numeric(model$times$final[1]))
}

# Coarse CV grid
get.all.results <- function(ds, fit.svmlight=TRUE, fit.small.grid=TRUE){
  subset <- c(1,2,4,5)
  
  C <- 10^(-6:7)
  C_poly <- 10^(-4:3)
  Gamma <- 10^(-8:8)
  degree <- c(2,3)
  
  if(fit.small.grid){
    C <- 10^(-5:4)
    C_poly <- 10^(-4:3)
    Gamma <- 10^(-6:7)
    degree <- c(2,3)
  }
  
  model.names=list("gmum.r::svm.radial", "gmum.r::svm.linear", "gmum.r::svm.poly", 
                   "kernLab::svm.radial", "kernLab::svm.linear", "kernLab::svm.poly")
  model.calls=list(caret.gmumSvmRadial, caret.gmumSvmLinear, caret.gmumSvmPoly,
                   "svmRadial", "svmLinear", "svmPoly")
  model.tuneGrids = list(expand.grid(C=C, gamma=Gamma), expand.grid(C=C), 
                         expand.grid(C=C_poly, gamma=Gamma, degree=degree, coef0=c(0)),
                         expand.grid(C=C, sigma=Gamma), expand.grid(C=C), 
                         expand.grid(C=C_poly, scale=Gamma,  degree=degree))
  
  
  
  model.names=model.names[subset]
  model.calls=model.calls[subset]
  model.tuneGrids = model.tuneGrids[subset]
  
  df <- as.data.frame(as.matrix(ds$x))
  df$y = ds$y
  df$name
  
  model.results <- list()
  for(model.index in 1:length(model.names)){
    name <- model.names[[model.index]]
    if(grepl("^gmum", name)){
      print("Fitting")
      name.experiment <- paste("libsvm_",name,ds$name, sep="")
      print(name.experiment)
      model.results[[name.experiment]] = 
             get.results(name.experiment, df, "y", model.tuneGrids[[model.index]], model.calls[[model.index]], verbosity=0, lib='libsvm')
      if(fit.svmlight){
        print("Fitting")
        name.experiment <- paste("svmlight_",name,ds$name, sep="")
        print(name.experiment)
        model.results[[name.experiment]] = 
          get.results(name.experiment, df, "y", model.tuneGrids[[model.index]], model.calls[[model.index]], verbosity=0, lib='svmlight')
      }
    }else{
      print("Fitting")
      print(name)
      print(ds$name)
      model.results[[name]] = 
        get.results(name, df, "y", model.tuneGrids[[model.index]], model.calls[[model.index]])
    }
    save(list="model.results", file=paste("comparative.tests.", ds$name, model.index, ".RData", sep=""))
  }

  M <- data.frame(trainAcc=sapply(model.results, function(x) x$trainAcc),
                  trainAccStd=sapply(model.results, function(x) x$trainAccStd),
                  testAcc=sapply(model.results, function(x) x$testAcc),
                  trainTime=sapply(model.results, function(x) x$trainTime),
                  finalModelTrainTime=sapply(model.results, function(x) x$finalModelTime)
  )
  row.names(M) <- sapply(model.results, function(x) x$name)
  
  save(list="M", file=paste("comparative.tests.", ds$name, ".RData", sep=""))
  
  M
}

for(ds in datasets){
  M <- get.all.results(ds, fit.svmlight = TRUE, fit.small.grid = TRUE)
}





