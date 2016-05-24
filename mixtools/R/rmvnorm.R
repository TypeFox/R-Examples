# Alternative version of rmvnorm to eliminate dependence of mixtools
# on additional package 'mvtnorm'

# Uses eigen decomposition assuming symmetric==TRUE.  Don't know how efficient
# this might be relative to other approaches, but some suggest this is among
# the most efficient methods to find a matrix square root.

rmvnorm <- function(n, mu=NULL, sigma=NULL) {
  if (is.null(mu)) {
    if (is.null(sigma)) {
      return(rnorm(n))
    } else {
      mu = rep(0, nrow(sigma))
    }
  } else if (is.null(sigma)) {
    sigma=diag(length(mu))
  }
  lmu <- length(mu)
  if (lmu != nrow(sigma) || lmu != ncol(sigma)) 
    stop("length of mu must equal nrow and ncol of sigma")
  e <- eigen(sigma, symmetric=TRUE)
  if (any(e$val<0)) 
    stop("Numerically negative definite covariance matrix")
  t(mu + e$vec %*% (t(e$vec) * sqrt(e$val)) %*% matrix(rnorm(n*lmu), lmu, n))  
}



