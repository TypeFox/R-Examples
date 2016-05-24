################################
#### Random values generation from a multivariate normal 
#### distribution on the simplex
#### Tsagris Michail 02/2016
#### mtsagris@yahoo.gr
################################

rcompnorm <- function(n, m, s, type = "alr") {
  ## n is the sample size
  ## m is the mean vector in R^d
  ## s is the covariance matrix in R^d
  ## type is either alr or ilr
  x <- rmvnorm(n, m, s)
  if (type == "alr") {
    y <- cbind( 1, exp(x) )
  } else {
    D <- ncol(x)
    y <- x %*% helm(D + 1)
    y < exp(y)
  } 
  y / rowSums(y)
}
