################################
#### Discrminant analysis for directional data
#### assuming a von Mises-Fisher distribution
#### Tsagris Michail 03/2014
#### mtsagris@yahoo.gr
#### References: J. E. Morris and P. J. Laycock (1974)
#### Discriminant Analysis of Directional Data (Biometrika)
################################

vmfda.pred <- function(xnew, x, ina) {
  ## xnew is the new observation(s)
  ## x is the data set
  ## ina is the group indicator variable
  x <- as.matrix(x)
  x <- x / sqrt( rowSums(x ^ 2) )
  xnew <- as.matrix(xnew)
  if (ncol(xnew) == 1)  xnew <- t(xnew)
  xnew <- xnew / sqrt( rowSums(xnew ^ 2) )
  p <- ncol(x)  ## dimensonality of the data
  ina <- as.numeric(ina)
  g <- max(ina)
  mesi <- matrix(nrow = g, ncol = p)
  k <- numeric(g)
  nu <- nrow(xnew)
  mat <- matrix(nrow = nu, ncol = g)
  est <- numeric(nu)
  for (j in 1:g) {
    da <- vmf(x[ina == j, ])  ## estimates the parameters of the vMF
    mesi[j, ] <- da$mu  ## mean direction
    k[j] <- da$k  ## concentration
  }
  for (j in 1:g) {
    mat[, j] <- (p/2 - 1) * log(k[j]) + k[j] * xnew %*% mesi[j, ] - 0.5 * p *
    log(2 * pi) - log(besselI(k[j], p/2 - 1, expon.scaled = TRUE)) - k[j]
  }
  apply(mat, 1, which.max)
}
