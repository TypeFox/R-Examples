################################
#### Random values generation from a multivariate t distribution
#### Tsagris Michail 07/2014
#### mtsagris@yahoo.gr
################################

rmvt <- function(n, mu, sigma, v) {
  ## n is the sample size
  ## mu is the mean vector
  ## sigma is the covariance matrix
  ## sigma does not have to be of full rank
  ## v is the degrees of freedom
  p <- length(mu)
  x <- matrix(rnorm(n * p), ncol = p)
  w <- sqrt(v / rchisq(n, v))
  eig <- eigen(sigma)
  lam <- eig$values
  vec <- eig$vectors
  B <- vec %*% diag( sqrt(lam) ) %*% t(vec)
  w * tcrossprod(x, B) + rep( mu, rep(n, p) )
}