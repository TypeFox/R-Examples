# --------------------------------------------------------------------
# Method/Function: predictClassProb
# --------------------------------------------------------------------

## Description:
##  - predicted probabilities for newdata based on a fitted model  
##    if a predict-mehtod exists this is just a wrapper.
##
## Arguments:
##   obejct: a fitted model
##   newdata: test data set 

## RETURNS: matrix of predicted probabilities (row=subjects, col=classes)

predictClassProb <- function(object,...){
  UseMethod("predictClassProb",object)
}

predictClassProb.randomForest <- function(object, newdata){
  stopifnot(!missing(newdata))
  P <- predict(object,newdata=newdata,type="prob")
  P
}

predictClassProb.mvr <- function(object, newdata, ncomp = 5){
	stopifnot(!missing(newdata))
	P <- predict(object,ncomp=ncomp,newdata=newdata,type="response")
	P[,,1]
}
 
predictClassProb.autoElasticNet <- function(object, newdata,...){
  stopifnot(!missing(newdata))
  P <- predict(object$enet, newx=newdata, type="response", s=object$Lambda)
  P[,,1]
}
