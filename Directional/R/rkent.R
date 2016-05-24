################################
#### Simulating from a Kent distribution
#### Tsagris Michail 04/2016
#### mtsagris@yahoo.gr
#### References: A new method to simulate the Bingham and related distributions in
#### directional data analysis with applications
#### Kent J.T., Ganeiber A.M. and Mardia K.V. (2013)
#### http://arxiv.org/pdf/1310.8110v1.pdf
################################

rkent <- function(n, k, m, b) {
  ## n is the required sample size
  ## k is the concentration parameter, the Fisher part
  ## m is the mean vector, the Fisher part
  ## b is the ovalness parameter

  m0 <- rnorm(3)
  m0 <- m0 / sqrt( sum(m0^2) )
  m <- m / sqrt( sum(m^2) )
  a <- rotation(m0, m)
  A <- diag(c(-b, 0, b))
  x <- rfb(n, k, m0, A)  ## simulated values with mean direction equal to m
  x %*% t(a)  ## simulated values with mean direction equal to m

}
