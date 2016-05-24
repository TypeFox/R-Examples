glinternet = function(X, Y, numLevels, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, screenLimit=NULL, numToFind=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000, verbose=FALSE, numCores=1){

                                        #get call and family
  thisCall = match.call()
  family = match.arg(family)
  
                                        #make sure inputs are valid
  n = length(Y)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  stopifnot(n==nrow(X), pCat+pCont==ncol(X), family=="gaussian"||family=="binomial")
  if (family=="binomial" && !all(Y %in% 0:1)) stop("Error:family=binomial but Y not in {0,1}")

                                        #separate into categorical and continuous parts
  if (pCont > 0) Z = as.matrix(apply(as.matrix(X[, numLevels == 1]), 2, standardize))
  else Z = NULL
  if (pCat > 0){
    catIndices = which(numLevels > 1)
    levels = numLevels[catIndices]
    Xcat = as.matrix(X[, catIndices])
  }
  else {
    levels = NULL
    Xcat = NULL
  }
  
                                        #compute variable norms
  res = Y - mean(Y)
  candidates = get_candidates(Xcat, Z, res, n, pCat, pCont, levels, screenLimit, numCores=numCores)
 
 
                                        #lambda grid if not user provided
  if (is.null(lambda)) lambda = get_lambda_grid(candidates, nLambda, lambdaMinRatio)
  else {
    stopifnot(min(lambda) > 0)
    if (any(diff(lambda) > 0)) stop("Error: input lambda sequence is not monotone decreasing.")
    lambdaMax = max(get_lambda_grid(candidates, nLambda, lambdaMinRatio))
    nLambda = length(lambda)
    if (nLambda == 1){
      lambda = sort(c(lambda, lambdaMax), decreasing=TRUE)
      nLambda = 2
    } else if (lambda[1] > lambdaMax) {
        cat("Info: some lambda values are larger than necessary, modifying ...\n")
        idx = max(which(lambda >= lambdaMax))
        if (idx == nLambda) {
            stop("All lambda values result in zero estimates.")
        }
        lambda = c(lambdaMax, lambda[(idx+1):nLambda])
        nLambda = length(lambda)
    }
  }
  
                                        #initialize storage for results
  fitted = matrix(mean(Y), n, nLambda)
  activeSet = vector("list", nLambda)
  betahat = vector("list", nLambda)
  betahat[[1]] = ifelse(family=="gaussian", mean(Y), -log(1/mean(Y)-1))
  objValue = rep(0, nLambda)
  objValue[1] = ifelse(family=="gaussian", sum(res^2)/(2*n), -mean(Y)*betahat[[1]]+log(1/(1-mean(Y))))

                                        #ever-active set + sequential strong rules + group lasso
  for (i in 2:nLambda){
    if (verbose) cat("lambda ", i, ": ", lambda[i], "\n")
    activeSet[[i]] = strong_rules(candidates, lambda[i], lambda[i-1])
    betahat[[i]] = initialize_betahat(activeSet[[i]], activeSet[[i-1]], betahat[[i-1]], levels)
    while (TRUE){
      #group lasso on strong set
      solution = group_lasso(Xcat, Z, Y, activeSet[[i]], betahat[[i]], levels, lambda[i], family, tol, maxIter)
      activeSet[[i]] = solution$activeSet
      betahat[[i]] = solution$betahat
      res = solution$res
      objValue[i] = solution$objValue
      #check kkt conditions on the rest
      check = check_kkt(Xcat, Z, res, n, pCat, pCont, levels, candidates, activeSet[[i]], lambda[i], numCores)
      candidates$norms = check$norms
      if (check$flag) break
      betahat[[i]] = initialize_betahat(check$activeSet, activeSet[[i]], betahat[[i]], levels)
      activeSet[[i]] = check$activeSet
    }
    #update the candidate set if necessary
    if (!is.null(screenLimit) && (screenLimit<pCat+pCont) && i<nLambda) candidates = get_candidates(Xcat, Z, res, n, pCat, pCont, levels, screenLimit, activeSet[[i]], candidates$norms, numCores)
    #get fitted values
    fitted[, i] = Y - res
    #compute total number of interactions found
    if (!is.null(numToFind)){
      numFound = sum(sapply(activeSet[[i]][3:5], function(x) ifelse(is.null(x), 0, nrow(x))))
      if (numFound >= numToFind) break
    }
  }

  #rescale betahat
  Z = as.matrix(X[, numLevels==1])
  betahatRescaled = lapply(1:i, function(j) rescale_betahat(activeSet[[j]], betahat[[j]], Xcat, Z, levels, n))

  output = list(call=thisCall, fitted=fitted[, 1:i], lambda=lambda[1:i], objValue=objValue, activeSet=activeSet[1:i], betahat=betahatRescaled[1:i], numLevels=numLevels, family=family)
  class(output) = "glinternet"
  return (output)
}
