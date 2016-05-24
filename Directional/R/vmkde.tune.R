
################################
#### Kernel density estimation of circular data with a von Mises kernel
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### Tuning the bandwidth
################################

vmkde.tune <- function(u, low = 0.1, up = 1, rads = TRUE) {
  ## u is the data
  n <- length(u)  ## sample size
  ## if the data are in degrees we transform them into radians
  if (rads == FALSE)  u <- u/180 * pi
  disa <- as.matrix(dist(u))
  funa <- function(h) {
  A <- exp( cos(disa) / h^2 )
    diag(A) <- NA  ## we do not want the diagonal elements
    f <- rowSums( A, na.rm = TRUE )/( (n - 1) * 2 * pi * besselI(1/h^2, 0) )
    mean( log(f) )
    }
  bar <- optimize(funa, c(low, up), maximum = TRUE)
  res <- c( bar$maximum, bar$objective )
  names(res) <- c("Optimal h", "cv")
  res
}
