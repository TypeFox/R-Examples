# Alternative version of dmvnorm to eliminate dependence of mixtools
# on additional package 'mvtnorm'
# Written (hopefully) to be more efficient than mvtnorm version, which uses both
# a call to "eigen" and a call to "mahalanobis", by using only a single
# call to the more efficient "qr" (read "Note" under ?qr)

# Note:  These functions assume that each ROW of y is a separate position vector.
# i.e., y is assumed to be nxd, where d=dimension
dmvnorm <- function(y, mu=NULL, sigma=NULL) {
  exp(logdmvnorm(y, mu=mu, sigma=sigma))
}


logdmvnorm <- function(y, mu=NULL, sigma=NULL) {
  if (is.vector(y)) 
    y <- matrix(y, nrow=1)
  d <- ncol(y)
  if (!is.null(mu))
    y <- sweep(y, 2, mu, '-')
  if (is.null(sigma))
    sigma <- diag(d)
  k <- d * 1.8378770664093454836 # that constant is log(2*pi)
  a <- qr(sigma)
  logdet <- sum(log(abs(diag(a$qr))))
  if(nrow(y)==1) 
    mahaldist <- as.vector(y %*% qr.solve(a,t(y)))
  else
    mahaldist <- rowSums((y %*% qr.solve(a)) * y)
  -0.5*(mahaldist + logdet + k)
}

