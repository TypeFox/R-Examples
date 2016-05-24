# Copyright Giovanni Petris <GPetris@uark.edu> 2001
#
###
### Fit a Poisson Cluster Process
###
pcp <- function(point.data, poly.data, h0=NULL, expo=0.25, n.int=20) {
  ## point.data: a points object
  ## poly.data: a polygon enclosing the study region
  ## h0: upper bound of integration in the criterion function
  ## expo: exponent in the criterion function
  ## n.int: number of intervals used to approximate the integral
  ##   in the criterion function with a sum 
  if (is.null(h0)) {
    dsq <- dsquare(point.data, point.data)
    h0 <- sqrt(max(dsq)/3)
  }
  h <- h0 / 20 * 1:20
  ## Compute K hat
  K.hat <- khat(point.data, poly.data, h)
  ## Define a function that computes K(h;theta)
  ## theta[1] = log(sigma^2),
  ## theta[2] = log(rho)
  K <- function(h, theta) {
    theta <- exp(theta)
    pi*h^2 + (1 - exp(-h^2/(4*theta[1])))/theta[2]
  }
  ## Define a function that evaluates the criterion
  D <- function(theta) {
    K.values <- K(h, theta)
    sum((K.hat^expo - K.values^expo)^2)
  }
  ## Minimize the criterion
  fit <- optim(c(0,0), D)
  fit$par <- exp(fit$par)
  names(fit$par) <- c("s2","rho")
  fit
}

