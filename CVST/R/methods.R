constructSVRLearner = function() {
  learn.svr = function(data, params) {
    #require(kernlab)
    stopifnot(isRegression(data))
    kpar=params[setdiff(names(params), c("kernel", "nu", "C"))]
    return(ksvm(data$x, data$y, kernel=params$kernel, kpar=kpar, type="nu-svr", nu=params$nu, C=params$C / getN(data), scale=FALSE))
  }
  
  predict.svr = function(model, newData) {
    stopifnot(isRegression(newData))
    return(predict(model, newData$x))
  }
  return(constructLearner(learn.svr, predict.svr))
}

constructSVMLearner = function() {
  learn.svm = function(data, params) {
    #require(kernlab)
    stopifnot(isClassification(data))    
    kpar=params[setdiff(names(params), c("kernel", "nu"))]
    return(ksvm(data$x, data$y, kernel=params$kernel, kpar=kpar, type="nu-svc", nu=params$nu, scale=FALSE))
  }
  
  predict.svm = function(model, newData) {
    stopifnot(isClassification(newData))    
    return(predict(model, newData$x))
  }
  return(constructLearner(learn.svm, predict.svm))
}

constructKlogRegLearner = function() {
  learn.klogreg = function(data, params) {
    #require(kernlab)
    stopifnot(isClassification(data))    
    # convert the factor to numeric 0/1
    if (nlevels(data$y) > 2) {
      stop("klogreg does not support multiclass experiments")
    }
    y = (data$y != levels(data$y)[1]) + 0
    kpar = params[setdiff(names(params), c("kernel", "lambda", "tol", "maxiter"))]
    kernel = do.call(params$kernel, kpar)
    model = .klogreg(data$x, kernel, y, getN(data) * params$lambda, params$tol, params$maxiter)
    model$yLevels = levels(data$y)
    return(model)
  }
  
  predict.klogreg = function(model, newData) {
    stopifnot(isClassification(newData))    
    pred = .klogreg.predict(model, newData$x)
    f = factor(pred, c("0", "1"), model$yLevels, ordered=FALSE)
    return(f)
  }
  return(constructLearner(learn.klogreg, predict.klogreg))
}

constructKRRLearner = function() {
  learn.krr = function(data, params) {
    #require(kernlab)
    stopifnot(isRegression(data))
    kpar = params[setdiff(names(params), c("kernel", "lambda"))]
    kernel = do.call(params$kernel, kpar)
    return(.krr(data$x, kernel, data$y, getN(data) * params$lambda))
  }
  
  predict.krr = function(model, newData) {
    stopifnot(isRegression(newData))
    return(as.matrix(.krr.predict(newData$x, model)))
  }
  return(constructLearner(learn.krr, predict.krr))
}

.krr = function(data, kernel, y, lambda) {
  #require(kernlab)
  #require(Matrix)
  K = kernelMatrix(kernel, data)
  N = nrow(K)
  alpha = solve(Matrix(K + diag(lambda, N))) %*% y
  return(list(data=data, kernel=kernel, alpha=alpha))
}

.krr.predict = function(newData, krr) {
  #require(kernlab)
  k = kernelMatrix(krr$kernel, newData, krr$data)
  return(k %*% krr$alpha)
}

.klogreg = function(data, kernel, labels, lambda, tol, maxiter) {
  # labels should be 0/1
  #require(kernlab)
  #require(Matrix)
  K = Matrix(kernelMatrix(kernel, data)@.Data)
  N = nrow(K)
  alpha = rep(1/N, N)
  iter = 1
  while (TRUE) {
    Kalpha = as.vector(K %*% alpha)
    spec = 1 + exp(-Kalpha)
    pi = 1 / spec
    diagW = pi * (1 - pi)
    e = (labels - pi) / diagW
    q = Kalpha + e
    theSol = try(solve(K + lambda * Diagonal(x=1/diagW), q))
    if (class(theSol) == "try-error") {
      break
    }
    alphan = as.vector(theSol)
    if (any(is.nan(alphan)) || all(abs(alphan - alpha) <= tol)) {
      break
    }
    else if (iter > maxiter) {
      cat("klogreg:maxiter!")
      break
    }
    else {
      alpha = alphan
      iter = iter + 1
    }
  }
  return(list(data=data, kernel=kernel, alpha=as.vector(alpha), pi=pi))
}

.klogreg.predict = function(klogreg, newData) {
  #require(kernlab)
  K = kernelMult(klogreg$kernel, newData, klogreg$data, klogreg$alpha)
  pi = 1 / (1 + exp(-as.vector(K)))
  return((pi >= .5) + 0)
}
