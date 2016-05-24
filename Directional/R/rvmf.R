################################
#### Simulating from a von Mises-Fisher distribution 
#### Tsagris Michail 10/2013 
#### mtsagris@yahoo.gr
#### References: Wood ATA (1994)
#### Simulation of the von Mises Fisher distribution (Communications in Statistics-Simulation) and
#### Inderjit S. Dhillon and Suvrit Sra (2003)
#### Modeling Data using Directional Distributions (Technical report, The University of Texas at Austin)
################################

rvmf <- function(n, mu, k) {
  ## n is the sample size
  ## mu is the mean direction and
  ## k is the concentration parameter
  ## n is the sample size
  d <- length(mu)  ## the dimensions
  if (k > 0) {
  mu <- mu/sqrt(sum(mu^2))  ## the mean direction
  ini <- c(numeric(d - 1), 1)  ## mean direction is set to (0, ..., 0, 1)
  b <- (-2 * k + sqrt(4 * k^2 + (d - 1)^2))/(d - 1)
  x0 <- (1 - b)/(1 + b)
  S <- matrix(nrow = n, ncol = d)
  m <- 0.5 * (d - 1)
  c <- k * x0 + (d - 1) * log(1 - x0^2)
  for (i in 1:n) {
    t <- -1000
    u <- 1
    while (t - c < log(u)) {
      z <- rbeta(1, m, m)
      u <- runif(1)
      w <- (1 - (1 + b) * z)/(1 - (1 - b) * z)
      t <- k * w + (d - 1) * log(1 - x0 * w)
    }
    v1 <- rnorm(d - 1)   
    v <- v1 / sqrt( sum(v1^2) )
    S[i, ] <- c(sqrt(1 - w^2) * v, w)
  }
  A <- rotation(ini, mu)  ## calculate the rotation matrix
  ## in order to rotate the initial mean direction from ini to mu
  x <- tcrossprod(S, A)  ## the x has direction mu 
  } else {  ## uniform distribution
  ## requires MASS if k = 0
  x1 <- MASS::mvrnorm( n, numeric(d), diag(d) )
  x <- x1 / sqrt( rowSums(x1^2) )
  } 
  x
}