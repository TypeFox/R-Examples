coef.logitchoice = function(object, lambdaIndex=NULL, ...){
  
  if (is.null(lambdaIndex)) {
    lambdaIndex = 1:length(object$lambda)
  }
  stopifnot(min(lambdaIndex)>0 && max(lambdaIndex)<=length(object$lambda))
  
  as.matrix(object$betahat[, lambdaIndex])
}