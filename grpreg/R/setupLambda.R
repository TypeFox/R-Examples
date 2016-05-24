setupLambda <- function(X, y, group, family, penalty, alpha, lambda.min, nlambda, group.multiplier) {
  ## Fit to unpenalized covariates
  n <- length(y)
  K <- table(group)
  K1 <- if (min(group)==0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  if (K1[1]!=0) {
    fit <- glm(y~X[, group==0], family=family)
  } else fit <- glm(y~1, family=family)

  ## Determine lambda.max
  if (family=="gaussian") {
    r <- fit$residuals
  } else {
    w <- fit$weights
    if (max(w) < 1e-4) stop("Unpenalized portion of model is already saturated; exiting...")
    r <- residuals(fit, "working")*w
  }
  if (strtrim(penalty,2)=="gr") {
    zmax <- .Call("maxgrad", X, r, K1, as.double(group.multiplier)) / n
  } else {
    zmax <- .Call("maxprod", X, r, K1, as.double(group.multiplier)) / n
  }
  lambda.max <- zmax/alpha
  
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
  } else {
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  }
  lambda
}

setupLambda.gBridge <- function(X, y, group, family, alpha, lambda.min, lambda.max, nlambda, gamma, group.multiplier)
{
  ## Fit to unpenalized covariates
  n <- length(y)
  ind <- which(group!=0)
  if (length(ind)!=length(group)) {
    fit <- glm(y~X[, group==0], family=family)
  } else {
    fit <- glm(y~1, family=family)
  }
  
  ## Guess lambda.max
  if (missing(lambda.max)) {
    if (family=="gaussian") {
      z <- crossprod(X[,ind], fit$residuals) / n
      a <- .35
    } else {
      z <- crossprod(X[,ind], fit$weights * residuals(fit, "working")) / n
      a <- .2
    }
    maxGradient <- tapply(abs(z), group[ind],max)*a^(1-gamma)/gamma
    lambda.max <- max(maxGradient/group.multiplier) / alpha
  }
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), len=nlambda-1)),0)                  
  } else {
    lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  }
  return(rev(lambda))
}
