logitchoice = function(X, Y, grouping, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, tol=1e-3, alpha=0.8, maxIter=5000, verbose=FALSE, numCores=1) {
    thisCall = match.call()

    numGroups = length(unique(grouping))
    groupSizes = sapply(unique(grouping), function(x) sum(grouping == x))
   
    #check inputs
    n = length(Y)
    stopifnot(n == sum(groupSizes))
    stopifnot(n == nrow(X))
    stopifnot(all(unique(Y) %in% c(0,1)))
    stopifnot(sum(Y) == length(groupSizes))
    stopifnot(tol > 0)
    stopifnot(maxIter > 0)
    stopifnot(numCores > 0)
    p = ncol(X)

    #get lambda values if not already supplied
    if (is.null(lambda)) {
      lambdaMax = 10
      lambdaMin = lambdaMinRatio * lambdaMax
      f = seq(0,1,1/(nLambda-1))
      lambda = lambdaMax^(1-f) * lambdaMin^f
    }
    else {
      stopifnot(all(lambda > 0))
      stopifnot(all(diff(lambda) <= 0))
    }
    nLambda = length(lambda)
    
    betahatMatrix = matrix(0, p , nLambda)
    probabilityMatrix = matrix(0, n, nLambda)
    residualMatrix = matrix(0, n, nLambda)
    numIters = rep(0, nLambda)
    objValues = list()
    betahat = rep(0, p)
    Yhat = rep(1/groupSizes, groupSizes)
    res = Y - Yhat
    for (i in 1:nLambda) {
      if (verbose) {
        cat("lambda ", i, ": ", lambda[i], "\n")
      }
      fit = .Call("R_solver", X, Y, res, Yhat, betahat, lambda[i], groupSizes, n, p, numGroups, alpha, 0, 0, maxIter, tol, double(maxIter), numCores)
      betahat = fit$betahat
      Yhat = fit$yhat
      res = fit$res
      betahatMatrix[, i] = betahat
      probabilityMatrix[, i] = Yhat
      residualMatrix[, i] = res
      numIters[i] = fit$numIters
      objValues = c(objValues, list(fit$objValue[1:(fit$numIters+1)]))
    }
    
    output = list(call=thisCall, betahat=betahatMatrix, yhat=probabilityMatrix, residual=residualMatrix, lambda=lambda, objValues=objValues, numIters=numIters)
    class(output) = "logitchoice"
    return (output)
}