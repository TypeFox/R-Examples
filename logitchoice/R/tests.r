#rm(list = ls())
#dyn.load("logitchoice.so")

test_compute_gradient = function(p, numCores) {
    groupSizes = sample(6:18, 1000, replace=TRUE)
    numGroups = length(groupSizes)
    n = sum(groupSizes)
    X = matrix(rnorm(n*p), nrow=n)
    res = rnorm(n)
    beta = rnorm(p)
    lambda = 0.1
    baseline = -t(X) %*% res / numGroups + lambda * beta
    cat("start computing\n")
    cat(paste(paste(dim(X), collapse=" "), "\n", sep=""))
    test = .Call("R_compute_gradient", X, res, beta, lambda, n, p, numGroups, numCores, double(p))
    max(abs(baseline-test))
}

test_compute_probabilities = function(p) {
  groupSizes = sample(6:18, 1000, replace=TRUE)
  numGroups = length(groupSizes)
  n = sum(groupSizes)
  X = matrix(rnorm(n*p), nrow=n)
  beta = rnorm(p)
  indices = c(0, cumsum(groupSizes))
  baseline = exp(X %*% beta)
  for (i in 1:numGroups) {
    range = (indices[i]+1) : indices[i+1]
    baseline[range] = baseline[range] / sum(baseline[range])
  }
  cat("start computing\n")
  cat(paste(paste(dim(X), collapse=" "), "\n", sep=""))
  test = .Call("R_compute_probabilities", X, beta, groupSizes, n, p, numGroups, double(n))
  max(abs(baseline-test))
}

test_solver = function(p, lambda=0, alpha=0.8, maxIter=5000, tol=1e-5, numCores=1) {
  groupSizes = sample(6:18, 1000, replace=TRUE)
  numGroups = length(groupSizes)
  n = sum(groupSizes)
  X = matrix(rnorm(n*p), nrow=n)
  X = scale(X)
  Y = rep(0, n)
  Y[cumsum(groupSizes)] = 1
  beta = rep(0, p)
  probabilities = unlist(sapply(groupSizes, function(x) rep(1/x, x)))
  res = Y - probabilities
  fit = .Call("R_solver", X, Y, res, probabilities, beta, lambda, groupSizes, n, p, numGroups, alpha, 0, 0, maxIter, tol, double(maxIter), numCores)
  c(max(abs(-t(X)%*%fit$res/numGroups + lambda*fit$betahat)), fit)
}

test_logitchoice = function(p) {
  groupSizes = sample(6:18, 1000, replace=TRUE)
  numGroups = length(groupSizes)
  n = sum(groupSizes)
  X = matrix(rnorm(n*p), nrow=n)
  X = scale(X)
  Y = rep(0, n)
  Y[cumsum(groupSizes)] = 1
  grouping = rep(1:numGroups, groupSizes)
  fit = logitchoice(X, Y, grouping)
  stopifnot(sum(fit$yhat) == numGroups*length(fit$lambda))
  list(yhat=max(abs(fit$yhat - predict.logitchoice(fit, X, grouping))), res=max(abs(fit$residual + fit$yhat - matrix(rep(Y, 50), n))))
}