coef.glinternet.cv = function(object, lambdaType=c("lambdaHat", "lambdaHat1Std"), ...){

  lambdaType = match.arg(lambdaType)
  
  bestIndex1Std = which(object$cvErr <= min(object$cvErr)+object$cvErrStd[which.min(object$cvErr)])
  if (lambdaType=="lambdaHat" || length(bestIndex1Std)==length(object$lambda)) return (extract_effects(object$betahat[[1]], object$activeSet[[1]], object$numLevels))

  idx = bestIndex1Std[1]
  extract_effects(object$glinternetFit$betahat[[idx]], object$glinternetFit$activeSet[[idx]], object$numLevels)
}
