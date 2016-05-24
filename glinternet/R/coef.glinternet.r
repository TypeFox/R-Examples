coef.glinternet = function(object, lambdaIndex=NULL, ...){

  if (is.null(lambdaIndex)) lambdaIndex = 1:length(object$lambda)
  stopifnot(min(lambdaIndex)>0 && max(lambdaIndex)<=length(object$lambda))
  
  lapply(1:length(object$lambda), function(x) extract_effects(object$betahat[[x]], object$activeSet[[x]], object$numLevels))
}
