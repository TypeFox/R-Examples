predict.glinternet.cv = function(object, X, type=c("response", "link"), lambdaType=c("lambdaHat", "lambdaHat1Std"), ...){

  type = match.arg(type)
  lambdaType = match.arg(lambdaType)
  lambda = ifelse(lambdaType=="lambdaHat", object$lambdaHat, object$lambdaHat1Std)
  
  predict(object$glinternetFit, X, type, lambda)
}

  
