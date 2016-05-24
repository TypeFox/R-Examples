
setupLambda <- function(X, y, family, lambda.min.ratio, nlambda, penalty.factor) {
  n <- nrow(X)
  p <- ncol(X)
  
  ## Determine lambda.max
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    fit <- glm(y~X[, -ind], family=family)
  } else fit <- glm(y~1, family=family)
  if (family=="gaussian") {
    lambda.max  <- .Call("maxprod", X, fit$residuals, ind, penalty.factor) / n
  } else {
    lambda.max  <- .Call("maxprod", X, residuals(fit, "working") * fit$weights, ind, penalty.factor) / n
  }
  
  if (lambda.min.ratio==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.max),log(lambda.min.ratio*lambda.max),len=nlambda))
  lambda
}