lmm.simu <- function(tau, sigma2, K, eigenK = eigen(K), X, beta) {
  if(any(eigenK$values < -1e-5))
    stop("K is not positive")
  eigenK$values[eigenK$values < 0] <- 0
  # partie fixe
  xbeta <- if(!missing(X)) {
    if(nrow(X) != nrow(eigenK$vectors)) stop("Dimensions mismatch (K and X)")
    if(length(beta) != ncol(X)) stop("Dimensions mismatch (X and beta)")
    X %*% beta
  } else 0
  # partie aléatoire
  omega <- eigenK$vectors %*% rnorm( nrow(eigenK$vectors), sd = sqrt(tau*eigenK$values) ) 
  y <- xbeta + omega + rnorm( nrow(eigenK$vectors), sd=sqrt(sigma2) )
  return( list(y = y, omega = omega) )
}
