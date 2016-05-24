get_lambda_grid = function(candidates, nLambda, lambdaMinRatio){
  
  lambdaMax = max(sapply(candidates$norms, function(x) ifelse(is.null(x), 0, max(x))))
  lambdaMin = lambdaMinRatio * lambdaMax
  f = seq(0,1,1/(nLambda-1))
  lambda = lambdaMax^(1-f) * lambdaMin^f
}

