predict.glinternet = function(object, X, type=c("response", "link"), lambda=NULL, ...){

  type = match.arg(type)
  
  helper = function(activeSet, betahat, numLevels, family){
    if (is.null(activeSet)){
      if (family=="gaussian" || type=="link") return(rep(betahat, n))
      return(rep(1/(1+exp(-betahat)), n))
    }
    nVars = sapply(activeSet, function(x) if(is.null(x)) 0 else nrow(x))
    indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)
    linear = .Call("R_x_times_rescaled_beta", Xcat, Z, betahat, n, nVars, numLevels, indices$cat, indices$cont, indices$catcat, indices$contcont, indices$catcont, double(n))
    if (family=="gaussian" || type=="link") return(linear)
    return(1/(1+exp(-linear)))
  }
  
  stopifnot(type=="link" || type=="response")
  X = as.matrix(X)
  n = nrow(X)
  pCat = sum(object$numLevels > 1)
  pCont = length(object$numLevels) - pCat
  stopifnot(pCat+pCont==ncol(X))
  if (pCont > 0) Z = matrix(X[, object$numLevels==1], nrow=n) else Z = NULL
  if (pCat > 0){
    catIndices = which(object$numLevels > 1)
    Xcat = matrix(as.integer(X[, catIndices]), nrow=n)
    levels = object$numLevels[catIndices]
  }
  else {
    levels = NULL
    Xcat = NULL
  }

  #if lambda is null, predict on all the lambdas
  if (is.null(lambda)){  
    return(sapply(1:length(object$betahat), function(x) helper(object$activeSet[[x]], object$betahat[[x]], levels, object$family)))
  }
  #otherwise, match the lambda sequence with user's lambda
  idx = match(lambda, object$lambda, 0)
  if (any(idx==0)) stop("Input lambda sequence not used in model fitting.")
  return(sapply(idx, function(x) helper(object$activeSet[[x]], object$betahat[[x]], levels, object$family)))
}
