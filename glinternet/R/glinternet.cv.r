glinternet.cv = function(X, Y, numLevels, nFolds=10, lambda=NULL, nLambda=50, lambdaMinRatio=0.01, screenLimit=NULL, family=c("gaussian", "binomial"), tol=1e-5, maxIter=5000, verbose=FALSE, numCores=1){

                                        #get call and family
  thisCall = match.call()
  family = match.arg(family)

                                        #make sure inputs are valid
  n = length(Y)
  pCat = sum(numLevels > 1)
  pCont = length(numLevels) - pCat
  stopifnot(n==nrow(X), pCat+pCont==ncol(X), family=="gaussian"||family=="binomial")

                                        #compute lambda grid if not provided
  if (is.null(lambda)){
    if (pCont > 0) Z = apply(as.matrix(X[, numLevels==1]), 2, standardize) else Z = NULL
    if (pCat > 0){
      catIndices = which(numLevels > 1)
      Xcat = matrix(as.integer(X[, catIndices]), nrow=n)
      levels = numLevels[catIndices]
    }
    else {
      levels = NULL
      Xcat = NULL
    }
                                        #compute variable norms
    res = Y - mean(Y)
    candidates = get_candidates(Xcat, Z, res, n, pCat, pCont, levels, screenLimit, numCores=numCores)
                                        #lambda grid
    lambda = get_lambda_grid(candidates, nLambda, lambdaMinRatio)
  }
  else {
    nLambda = length(lambda)
  }
  
                                        #create the folds
  folds = sample(rep(1:nFolds, ceiling(n/nFolds)), n, replace=FALSE)
  
                                        #perform cv
  compute_loss = function(y, yhat, family){
    if (family == "gaussian") return (sum((y-yhat)^2)/(2*length(y)))
    yhat = sapply(yhat, function(x) min(max(1e-15, x), 1-1e-15))
    -(t(y)%*%log(yhat) + t(1-y)%*%log(1-yhat))/length(y)
  }
  loss = matrix(0, nFolds, nLambda)
  for (fold in 1:nFolds){
    testIndex = which(folds == fold)
    trainIndex = which(folds != fold)
    Xtrain = as.matrix(X[trainIndex, ])
    Ytrain = Y[trainIndex]
    Xtest = as.matrix(X[testIndex, ])
    Ytest = Y[testIndex]
    fitted = glinternet(Xtrain, Ytrain, numLevels, lambda, nLambda, lambdaMinRatio, screenLimit, numToFind=NULL, family, tol, maxIter, verbose, numCores)
    YtestHat = predict(fitted, Xtest, "response")
    loss[fold, ] = apply(YtestHat, 2, function(yhat) compute_loss(Ytest, yhat, family))
  }

                                        #compute cv errors and get minimum
  cv = apply(loss, 2, mean)
  cvStd = apply(loss, 2, sd)
  bestIndex1Std = which(cv <= min(cv)+cvStd[which.min(cv)])
  bestIndex = which.min(cv)
  if (length(bestIndex1Std) == nLambda) lambdaHat1Std = lambdaHat = lambda[bestIndex]
  else {
    lambdaHat1Std = lambda[bestIndex1Std[1]]
    lambdaHat = lambda[bestIndex]
  }

                                        #perform fit with chosen lambda
  fitted = glinternet(X, Y, numLevels, lambda[1:bestIndex], nLambda, lambdaMinRatio, screenLimit, family=family, tol=tol, maxIter=maxIter, verbose=verbose, numCores=numCores)

  output = list(call=thisCall, glinternetFit=fitted, fitted=fitted$fitted[, bestIndex], activeSet=fitted$activeSet[bestIndex], betahat=fitted$betahat[bestIndex], lambda=lambda, lambdaHat=lambdaHat, lambdaHat1Std=lambdaHat1Std, cvErr=cv, cvErrStd=cvStd, family=family, numLevels=numLevels, nFolds=nFolds)
  class(output) = "glinternet.cv"
  return (output)
}
