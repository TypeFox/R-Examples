################################
#### MLE estimation of the Kent distribution on the sphere
#### Tsagris Michail 05/2014 
#### mtsagris@yahoo.gr
#### References: John Kent(1982)
#### JRSSB, 44(1): 71-80.
#### The Fisher-Bingham distribution on the sphere
################################

kent.mle <- function(x) {
  ## x is the data in Euclidean coordinates
  n <- nrow(x)  ## sample size
  xbar <- colMeans(x)  ## mean vector
  xbar <- xbar/sqrt(sum(xbar^2)) ## mean direction
  u <- c( acos(xbar[1]), ( atan(xbar[3] / xbar[2]) + pi * I(xbar[2]<0) )
  %% (2 * pi) )
  ## u is the mean vector to latitude and longitude
  theta <- u[1]
  phi <- u[2]
  H <- matrix(c(cos(theta), sin(theta) * cos(phi), 
       sin(theta) * sin(phi), -sin(theta), cos(theta) * cos(phi), 
       cos(theta) * sin(phi), 0, -sin(phi), cos(phi)), ncol = 3)
  S <- (1/n) * t(x) %*% x
  B <- t(H) %*% S %*% H
  psi <- 0.5 * atan(2 * B[2, 3]/(B[2, 2] - B[3, 3]))
  K <- matrix(c(1, 0, 0, 0, cos(psi), sin(psi), 0, 
       -sin(psi), cos(psi)), ncol = 3)
  G <- H %*% K  ## The G matrix Kent describes, the A in our notation
  r1 <- sqrt(sum(xbar^2))
  lam <- eigen(B[-1, -1])$values
  r2 <- lam[1] - lam[2]
  ## the next function will be used to estimate the kappa and beta
  mle <- function(para) {
    ## maximization w.r.t. to k and b
    k <- para[1]
    b <- para[2]
    gam <- c(0, k, 0)
    lam <- c(0, -b, b)
    ckb <- fb.saddle(gam, lam)[3]
    g <-  -( -n * ckb + k * sum(x %*% G[, 1]) + b * sum((x %*% G[, 2])^2) - 
    b * sum( (x %*% G[, 3])^2) )
    g
  }
  ini <- vmf(x)$k
  ini <- c(ini, ini/2.1)  ## initial values for kappa and beta
  qa <- optim(ini, mle)
  para <- qa$par
  k <- para[1]
  b <- para[2]  ## the estimated parameters
  gam <- c(0, k, 0)
  lam <- c(0, -b, b)
  ckb <- as.numeric( fb.saddle(gam, lam)[3] )
  ## the line below calculates the log-likelihood
  l <-  -n * ckb + k * sum(x %*% G[, 1]) + b * sum((x %*% G[, 2])^2) - 
        b * sum((x %*% G[, 3])^2)
  para <- c(k, b)
  names(para) <- c("kappa", "beta")
  colnames(G) <- c("mean", "major", "minor")
  list(G = G, para = para, logcon = ckb, loglik = l)
}