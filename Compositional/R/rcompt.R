################################
#### Random values generation from a multivariate 
#### t distribution on the simplex
#### Tsagris Michail 02/2016
#### mtsagris@yahoo.gr
################################

rcompt <- function(n, m, s, dof, type = "alr") {
  ## n is the sample size
  ## m is the mean vector in R^d
  ## s is the covariance matrix in R^d
  ## dof is the degrees of freedom
  ## type is either alr or ilr
  x <- rmvt(n, m, s, dof)
  if (type == "alr") {
    y <- cbind( 1, exp(x) )
  } else {
    D <- ncol(x)
    y <- x %*% helm(D + 1)
    y < exp(y)
  } 
  y / rowSums(y)
}
