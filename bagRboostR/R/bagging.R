#' @import randomForest 
library(randomForest)

#' @export
bagging <- function(formula,data,test,m=5,ntree=500,mtry=NULL,trace=T) {
  outcome.label <- outcomeLabel(formula)
  mtry <- ifelse(is.null(mtry),floor(sqrt(ncol(data))),mtry)
  C <- matrix(nrow=nrow(test),ncol=m)
  n <- nrow(data)
  for (i in 1:m) {
    t <- data[sample(n, n, replace=T), ]
    t[,outcome.label]<-droplevels(t[,outcome.label]) # in case any outcomes are not sampled
    fit <- randomForest(formula,data=t,ntree=ntree,do.trace=trace)      
    C[,i] <- as.character(predict(fit,test))
  }
  apply(C,1,bagPrediction)
}

bagPrediction <- function(sample) {
  max.classes <- maxClasses(sample)
  ifelse(length(max.classes)==1,max.classes,sample(max.classes,1))
}

maxVoteCount <- function(sample) {
  sum(sample==names(which.max(table(sample))))
}

maxClasses <- function(sample) {
  names(which(table(sample)==maxVoteCount(sample)))
}



