rescale_betahat = function(activeSet, betahat, X, Z, levels, n){

  if (is.null(activeSet)) return(betahat)
  
  nVars = sapply(activeSet, function(x) if (is.null(x)) 0 else nrow(x))
  betaLen = length(betahat)
  indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)

  .Call("R_rescale_beta", X, Z, n, betahat, betaLen, nVars, levels, indices$cat, indices$cont, indices$catcat, indices$contcont, indices$catcont, double(betaLen))
}

  
