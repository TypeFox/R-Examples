################################
#### Contour plots of the von Mises-Fisher distribution on the sphere
#### Tsagris Michail 06/2014 
#### mtsagris@yahoo.gr
################################

vmf.contour <- function(k) {
  ## k is the concentration parameter
  rho <- pi/2  ## radius of the circular disc
  x <- y <- seq(-rho, rho, by = 0.01)
  n <- length(x)
  mat1 <- matrix(rep(x^2, n), ncol = n)
  mat2 <- matrix(rep(y^2, n), ncol = n, byrow = T)
  z <- mat1 + mat2
  theta <- sqrt(z)
  ind <- ( theta < rho )  ## checks if x^2+y^2 < rho^2
  ind[ ind == FALSE ] <- NA
  xa <- 0.5 * log(k) + k * cos(theta) - 1.5 * log(2 * pi) - 
  log(besselI(k, 0.5, expon.scaled = T)) - k
  mat <- exp(xa) * ind   
  contour(x, y, mat)
}