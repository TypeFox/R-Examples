.rmvcauchy <- function(n,
                       mean=rep(0, nrow(sigma)),
                       sigma=diag(length(mean)),
                       method=c("eigen", "svd", "chol")) {
  # a multivariate cauchy-like distribution
  # this is wider tailed and lacks the Gaussian's tendency to cluster
  # PS.  I know this isn't an actual covariance
  # PPS.  This is a direct ripoff of mvtnorm:rmvnorm
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != nrow(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  sigma1 <- sigma
  dimnames(sigma1) <- NULL
  if (!isTRUE(all.equal(sigma1, t(sigma1)))) {
    warning("sigma is numerically not symmetric")
  }
  method <- match.arg(method)
  if (method == "eigen") {
    ev <- eigen(sigma, symmetric = TRUE)
    if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% 
      t(ev$vectors)
  }
  else if (method == "svd") {
    sigsvd <- svd(sigma)
    if (!all(sigsvd$d >= -sqrt(.Machine$double.eps) * abs(sigsvd$d[1]))) {
      warning("sigma is numerically not positive definite")
    }
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
  }
  else if (method == "chol") {
    retval <- chol(sigma, pivot = TRUE)
    o <- order(attr(retval, "pivot"))
    retval <- retval[, o]
  }
  retval <- matrix(rcauchy(n * ncol(sigma)), nrow = n) %*% retval
  retval <- sweep(retval, 2, mean, "+")
  colnames(retval) <- names(mean)
  retval
}
