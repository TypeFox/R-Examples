"ftrend" <-
function(object, ...)
{
  if(length(object$xlevels) == 0) {
    stop("No factors in model")
  }
  
  xname <- names(object$xlevels)[1]
  if (!identical(object$contrasts[[1]], "contr.treatment")) {
    stop(paste("Treatment contrasts must be used for variable", xname))
  }
  xlevels <- object$xlevels[[1]]
  nlevels <- length(xlevels)
  coefnames <- paste(xname, xlevels, sep="")
  ncoef <- length(coef(object))
  
  if (!all(coefnames %in% names(coef(object)))) {
    stop("The model must not have an intercept term")
  }
  
  index1 <- match(coefnames, names(coef(object)))
  index2 <- (1:ncoef)[-index1]
  m <- length(index1)
  ncov <- length(index2)

  ## Centre the covariates according to Greenland et al (weighted mean = 0)
  X0 <- model.matrix(object)
  if (!is.null(object$weights)) {
    mu <- -apply(X0, 2, weighted.mean, object$weights )[index2]
  }
  else {
    mu <- -apply(X0, 2, mean)[index2]
  }
  mu.full <- rep(0, ncoef)
  mu.full[index2] <- mu
  X <- sweep(X0, 2, mu.full, "+")

  ## Information matrix with centred covariates
  if (!is.null(object$weights)) {
    J <- crossprod(X, sweep(X, 1, object$weights, "*"))
  }
  else {
    ## linear models
    J <- crossprod(X,X)
  }

  J11 <- J[index1, index1]
  J12 <- J[index1, index2]
  J21 <- J[index2, index1]
  J22 <- J[index2, index2]

  ## Variance matrix
  V <- solve(J)
  V11 <- V[index1, index1]
  V12 <- V[index1, index2]
  V21 <- V[index2, index1]
  V22 <- V[index2, index2]

  cal.V <- function(mu) {
    one <- as.matrix(rep(1,m))
    mu <- as.matrix(mu)
    return(V11 - one %*% t(mu) %*% V21 - V12 %*% mu %*% t(one) +
           matrix(1, m, m) * (t(mu) %*% V22 %*% mu)[1,1])
  }

  f <- function(mu) {
    V.mu <- cal.V(mu)
    # lambda is current vector of floating variances
    D <- sum(diag(V.mu)/lambda) - c(determinant(V.mu)$modulus) +
      sum(log(lambda)) - m

    S <- (1/sum(diag(J11)) + t(mu) %*% solve(J22) %*% mu)[1,1]
    
    grad1 <- - t(1/lambda) %*% V12
    grad2 <- + sum(1/lambda) * t(mu) %*% V22
    grad3 <- - t(mu) %*% solve(J22)/S
    attr(D,"gradient") <- 2 * (grad1 + grad2 +  grad3)

    H1 <- V22 * sum(1/lambda)
    H2 <- -solve(J22)/S
    b <- solve(J22) %*% mu/S
    H3 <- 2 * b %*% t(b)
    attr(D, "hessian") <- 2 * (H1 + H2 + H3)
    return(D)
  }

  ## Initial value of lambda
  lambda <- diag(V[index1,index1])

  ## Do the minimization
  nlm.out <- nlm(f, rep(0,ncov), check.analyticals=TRUE, ...)
  mu2 <- nlm.out$estimate

  ## Calculate parameter values and their covariance matrix
  ## if the covariates are appropriately centred
  coef <- coef(object)[index1] - c(crossprod(mu + mu2, coef(object)[index2]))
  return(list(coef=coef, vcov=cal.V(mu2)))
}
