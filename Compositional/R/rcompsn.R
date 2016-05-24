################################
#### Random values generation from a multivariate skew 
#### normal distribution on the simplex
#### Tsagris Michail 02/2016
#### mtsagris@yahoo.gr
################################

rcompsn <- function(n, xi, Omega, alpha, dp = NULL, type = "alr") {
  ## n is the sample size
  ## m is the mean vector in R^d
  ## s is the covariance matrix in R^d
  ## type is either alr or ilr
  if ( is.null(dp) ) {
    x <- rmsn(n = n, xi = xi, Omega = Omega, alpha = alpha)
  } else  x <- rmsn(n = n, dp = dp)
  if (type == "alr") {
    y <- cbind( 1, exp(x) )
  } else {
    D <- ncol(x)
    y <- x %*% helm(D + 1)
    y < exp(y)
  } 
  y / rowSums(y)
}
