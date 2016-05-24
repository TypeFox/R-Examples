#' Permutation- based Variable Importance Measures
#'
#' \code{variableImportance} produces permutation- based variable importance measures (currently only for binary classification models from the package \code{randomForest} and only for the performance measure AUROC)
#'
#' @param object A model. Currently only binary classification models from the package \code{randomForest}.
#' @param xdata A data frame containing the predictors for the model. 
#' @param ydata A factor containing the response variable.
#' @param CV Cross-validation. How many times should the data be permuted and the decrease in performance be calculated? Afterwards the mean is taken. CV should be higher for very small samples to ensure stability.
#' @param measure Currently only Area Under the Receiver Operating Characteristic Curve (AUROC) is supported.
#' @param sort Logical. Should the results be sorted from high to low?
#' @details Currently only binary classification models from \code{randomForest} are supported. Also, currently only AUROC is supported. Definition of MeanDecreaseAUROC: for the entire ensemble the AUROC is recorded on the provided xdata. The same is subsequently done after permuting each variable (iteratively, for each variable separately). Then the latter is subtracted from the former. This is called the Decrease in AUROC. If we do this for multiple CV, it becomes the Mean Decrease in AUROC.
#' @examples
#' #Prepare data
#' data(iris)
#' iris <- iris[1:100,]
#' iris$Species <- as.factor(ifelse(factor(iris$Species)=="setosa",0,1))
#' #Estimate model
#' library(randomForest)
#' ind <- sample(nrow(iris),50)
#' rf <- randomForest(Species~., iris[ind,])
#' #Obtain variable importances
#' variableImportance(object=rf, xdata=iris[-ind,names(iris) != "Species"],
#' ydata=iris[-ind,]$Species) 
#' @seealso \code{\link{parDepPlot}}
#' @return  A data frame containing the variable names and the mean decrease in AUROC
#' @author Authors: Michel Ballings, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
variableImportance <- function(object=NULL,xdata=NULL,ydata=NULL, CV=3, measure="AUROC", sort=TRUE){
  
  #auc before permutation
  pred <- predict(object,newdata=xdata,type="prob")[,2]
  
  auc <- numeric()
  auc <- AUC::auc(roc(as.numeric(pred),as.factor(ydata)))
  
  store <- list()
  

  for (i in seq_len(CV)) {
      
      
      store[[i]] <- numeric(ncol(xdata))
      
      #auc after permutation
      for (ii in 1:ncol(xdata)){
        xdata.after <- xdata
        xdata.after[,ii] <- sample(xdata[,ii])
        
        pred <- predict(object,newdata=xdata.after,type="prob")[,2]
     
        
        auc.after <-  AUC::auc(roc(as.numeric(pred),as.factor(ydata)))
        store[[i]][ii] <- auc - auc.after
      }
  
  }
  
  
  store <- data.frame(VarName=colnames(xdata),MeanDecreaseAUROC=rowMeans(do.call(cbind,store)))
  
  cat(paste("Measure used:", measure, "\n"))
  
  if (sort==TRUE) {
    store[order(-store[,2]),]
  } else {
    store
  }
}