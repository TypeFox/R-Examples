

################################
#### Kernel density estimation of circular data with a von Mises kernel
#### Tsagris Michail 2/2015
#### mtsagris@yahoo.gr
#### References: C.C. Taylor (2008)
#### Automatic bandwidth selection for circular density estimation
#### Computational Statistics \& Data Analysis and
#### Garcia-Portugues E. (2013)
#### Exact risk improvement of bandwidth selectors for kernel
#### density estimation with directional data
#### Electronic Journal of Statistics
################################

vm.kde <- function(u, h = NULL, thumb = "none", rads = TRUE) {
  ## u is the data
  ## h is the bandwidth you want
  ## thumb is either 'none' (defualt), or 'tay' (Taylor, 2008) or
  ## 'rot' (Garcia-Portugues, 2013)
  ## if the data are in degrees we transform them into radians

  if (rads == FALSE)  u <- u/180 * pi
  n <- length(u)  ## sample size
  disa <- dist(u, diag = TRUE, upper = TRUE)
  disa <- as.matrix(disa)

  if ( is.null(h) ) {

    if (thumb == "tay") {
      k <- circ.summary(u, rads = TRUE, plot = FALSE)$kappa
      h <- ( (4 * pi^0.5 * besselI(k, 0)^2)/(3 * n * k^2 *
      besselI(2 * k, 2)) )^(1/5)
    } else if (thumb == "rot") {
      k <- circ.summary(u, rads = TRUE, plot = FALSE)$kappa
      h <- ( (4 * pi^0.5 * besselI(k, 0)^2) /( k * n * ( 2 * besselI(2 * k, 1) +
      3 * k * besselI(2 * k, 2)) ) )^(1/5)
    } else if (thumb == "none") {
      h <- as.numeric( vmkde.tune(u, low = 0.1, up = 1)[1] )
    }

  } else h <- h

  f <- rowSums( exp( cos(disa)/h^2 ) ) / ( n * 2 * pi * besselI(1/h^2, 0) )

  list( h = h, f = as.vector(f) )
}
