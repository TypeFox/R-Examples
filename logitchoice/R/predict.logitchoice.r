predict.logitchoice = function(object, X, grouping, lambda=NULL, ...) {
  
  grouping = as.factor(grouping)
  groupSizes = sapply(levels(grouping), function(x) sum(grouping == x))
  numGroups = length(groupSizes)
  n = nrow(X)
  stopifnot(n == sum(groupSizes))
  
  #if lambda is null, predict on all the lambdas
  if (!is.null(lambda)) {
    idx = match(lambda, object$lambda, 0)
    if (any(idx==0)) stop("Input lambda sequence not used in model fitting.")
  } else {
    idx = 1:length(object$lambda)
  }
  result = as.matrix(exp(X %*% object$betahat[, idx]))
  indices = c(0, cumsum(groupSizes))
  for (i in 1:numGroups) {
    temp = as.matrix(result[(indices[i]+1):indices[i+1], ])
    result[(indices[i]+1):indices[i+1], ] = scale(temp, center=FALSE, scale=apply(temp, 2, sum))
  }
  return (result)
}