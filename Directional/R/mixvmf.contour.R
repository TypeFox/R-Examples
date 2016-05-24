################################
#### Contour plots of a spherical model based clustering using mixtures 
#### of von Mises-Fisher distributions
#### Tsagris Michail 4/2015 
#### mtsagris@yahoo.gr
#### References: Kurt Hornik and  Bettina Grun (2014)
#### movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions
#### http://cran.r-project.org/web/packages/movMF/vignettes/movMF.pdf
################################

mixvmf.contour <- function(u, mod) {
  ## u contains the data in latitude and longitude the first column is 
  ## the latitude and the second column is the longitude
  ## mod is a mix.vmf object
  u <- as.matrix(u)  ## makes sure u is a matrix
  n <- nrow(u)  ## sample size
  n1 <- 100
  n2 <- 100  ## n1 and n2 specify the number of points taken at each axis
  x1 <- seq(min(u[, 1]) - 5, max(u[, 1]) + 5, length = n1)  ## latitude
  x2 <- seq(min(u[, 2]) - 5, max(u[, 2]) + 5, length = n2)  ## longitude
  mat <- matrix(nrow = n1, ncol = n2)
  mu <- mod$param[, 1:3]  ## mean directions
  k <- mod$param[, 4]  ## concentration parameters
  p <- mod$param[, 5]  ## mixing probabilities
  g <- length(p)  ## how many clusters
  lika <- con <- numeric(g)
  for (l in 1:g) {
  con[l] <- 0.5 * log(k[l]) - 1.5 * log(2 * pi) - ( log(besselI(k[l], 0.5, 
  expon.scaled = T)) + k[l] )
  }
  for (i in 1:n1) {
    for (j in 1:n2) {
      #x <- c( cos(x1[i]) * cos(x2[j]), cos(x1[i]) * sin(x2[j]), sin(x2[j]) )
      x <- euclid( c(x1[i], x2[j]) )
      lika <- con + k * ( x %*% t(mu) )
      can <- sum( p * exp(lika) )
      if (abs(can) < Inf) {
        mat[i, j] <- can 
      } else  mat[i, j] <- NA 
    }
  }
  contour(x1, x2, mat, nlevels = 8, col = 4, xlab = "Latitude", ylab = "Longitude")
  points(u[, 1], u[, 2])
}