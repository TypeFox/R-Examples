################################
#### Spherical-spherical regression
#### Tsagris Michail 11/2013
#### mtsagris@yahoo.gr
#### References: Chang Ted (1986)
#### Spherical egession. Annals of statistics, 14(3): 907-924
################################

spher.reg <- function(y, x, rads = FALSE) {
  ## x is the independent variable
  ## y is the dependent variable
  ## The first row of both matrices is the latitude
  ## and the second is the longitude

  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)  ## sample size

  if ( ncol(x) == 2  &  ncol(y) == 2 ) {
    if (rads == FALSE) {
      x <- pi * x / 180  ## from degrees to rads
      y <- pi * y / 180
    }  ## from degrees to rads
    ## the first row of both matrices is the latitude and the second is the longitude
    ## the next two rows transform the data to Euclidean coordinates
    X <- cbind(cos(x[, 1]) * cos(x[, 2]), cos(x[, 1]) * sin(x[, 2]), sin(x[, 1]))
    Y <- cbind(cos(y[, 1]) * cos(y[, 2]), cos(y[, 1]) * sin(y[, 2]), sin(y[, 1]))
  } else if ( ncol(x) == 3  &  ncol(y) == 3 ) {
    X <- x
    Y <- y
  }

  XY <- crossprod(X, Y) / n
  b <- svd(XY)  ## SVD of the XY matrix
  A <- b$v %*% t(b$u)
  if (det(A) < 0) {
    b$u[, 3] <-  - b$u[, 3]
    A <- tcrossprod(b$v, b$u )
  }

  est <- tcrossprod(X, A)
  list(A = A, fitted = est)
}
