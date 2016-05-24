MVTMLE <- function (X, nu = 1, location = TRUE, eps = 1e-06, maxiter = 100) 
{
  Xnames <- colnames(X)
  X <- as.matrix(X)
  if (location == TRUE) {
    p <- ncol(X)
    if (nu < 1) 
      stop("'nu' must be >= 1 when 'location=TRUE'")
    Y <- cbind(X, 1)
    res0 <- MVTMLE0(Y, nu = nu - 1, prewhitened = FALSE, 
                    eps = eps, maxiter = maxiter)
    G <- res0$S/res0$S[p + 1, p + 1]
    mu <- G[1:p, p + 1]
    Sigma <- G[1:p, 1:p] - tcrossprod(mu)
    names(mu) <- Xnames
    rownames(Sigma) <- colnames(Sigma) <- Xnames
  }
  if (location == FALSE) {
    p <- ncol(X)
    res0 <- MVTMLE0(X, nu = nu, eps = eps, maxiter = maxiter)
    mu <- rep(0, p)
    Sigma <- res0$S
    names(mu) <- Xnames
    rownames(Sigma) <- colnames(Sigma) <- Xnames
  }
  if (is.numeric(location)) {
    X <- t(t(X)-location)
    res0 <- MVTMLE0(X, nu = nu, eps = eps, maxiter = maxiter)
    mu <- location
    Sigma <- res0$S
    names(mu) <- Xnames
    rownames(Sigma) <- colnames(Sigma) <- Xnames
  }
  
  if (res0$nG > eps) warning("MVTMLE did not converge", call. = FALSE)
  
  return(list(mu = mu, Sigma = Sigma, iter=res0$iter))
}


