group_lasso = function(X, Z, Y, activeSet, betahat, numLevels, lambda, family, tol, maxIter){

  #get size of each group, number of groups in each category, and total number of groups
  groupSizes = get_group_sizes(activeSet, numLevels)
  numGroups = sapply(activeSet, function(x) if (is.null(x)) 0 else nrow(x))
  totalGroups = sum(numGroups)

  #if active set is empty, just return estimate of the intercept
  if (totalGroups == 0){
    res = Y - mean(Y)
    if (family == "gaussian") betahat = mean(Y)
    else betahat = -log(1/mean(Y)-1)
    return(list(betahat=betahat, activeSet=activeSet, res=res))
  }

  n = length(Y)
  indices = lapply(activeSet, function(x) if (!is.null(x)) c(t(x)) else NULL)
  
  #fit and get new betahat, res, objValue
  fit = .Call("R_gl_solver", X, Z, Y, n, betahat[1], betahat[-1], double(n), double(n), numLevels, numGroups, indices$cat, indices$cont, indices$catcat, indices$contcont, indices$catcont, lambda, tol, 0.1, maxIter, 0, 0, double(maxIter), ifelse(family=="gaussian", 0, 1))
  res = fit$res
  objValue = fit$objValue

  #get the nonzero parts of betahat and update activeSet
  idx = .Call("R_retrieve_beta", fit$coefficients, groupSizes, totalGroups, integer(totalGroups), integer(length(fit$coefficients)))
  beta = c(fit$mu, fit$coefficients[idx$betaIdx != 0])
  range = c(0, cumsum(numGroups))
  activeSet = lapply(1:5, function(i){
    if (numGroups[i] > 0){
      index = which(idx$idx[(range[i]+1):range[i+1]] != 0)
      if (length(index) > 0) matrix(activeSet[[i]][index, ], nrow=length(index))
      else NULL
    }
    else NULL
  })
  names(activeSet) = c("cat", "cont", "catcat", "contcont", "catcont")

  #output
  list(betahat=beta, activeSet=activeSet, res=res, objValue=objValue)
}
  



