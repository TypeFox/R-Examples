################################
#### Saddlepoint approximations of the Fisher-Bingham distributions 
#### Tsagris Michail 02/2014 
#### mtsagris@yahoo.gr
#### References: Kume Alfred and Wood Andrew T.A. (2005)
#### Saddlepoint approximations for the Bingham and Fisher-Bingham normalizing constants (Biometrika)
################################

fb.saddle <- function(gam, lam) {
  ## gam is the parameters of the Fisher part
  ## lam is the eigenvalues of the matrix of the Bingham part
  lam <- sort(lam)  ## sorts the eigenvalues of the Bingham part
  mina <- min(lam)
  if (mina <= 0) {
      aaa <- abs(mina) + 1
      lam <- lam + aaa ## makes all the lambdas positive and greater than zero
  }  
  p <- length(gam)  ## dimensionality of the distribution
  para <- c(gam, lam)  ## the parameters of the Fisher-Bingham
  saddle.equat <- function(ta, para) {
    ## saddlepoint equation
    p <- length(para)/2
    gam <- para[1:p]
    lam <- para[ -c(1:p) ]
    f <- sum( 0.5/(lam - ta) + 0.25 * ( gam^2/(lam - ta)^2 ) ) - 1
    f
  }
  low <- lam[1] - 0.25 * p - 0.5 * sqrt(0.25 * p^2 + p * max(gam)^2)  ## lower bound
  up <- lam[1] - 0.25 - 0.5 * sqrt(0.25 + min(gam)^2)  ## not the exact upper 
  ## bound but a bit higher
  ela <- uniroot(saddle.equat, c(low, up), para = para, tol = 1e-08)
  tau <- ela$root  ## tau which solves the saddlepoint equation
  ### below are the derivatives of the cumulant generating function
  kfb <- function(j, gam, lam, ta) {
    if (j == 1) {
      kd <- sum( 0.5/(lam - ta) + 0.25 * ( gam^2/(lam - ta)^2 ) )
    } else if (j > 1) {
      kd <- sum( 0.5 * factorial(j - 1)/(lam - ta)^j + 0.25 * factorial(j) * 
        gam^2/(lam - ta)^(j + 1) )
    }
    kd
  }
  rho3 <- kfb(3, gam, lam, tau)/kfb(2, gam, lam, tau)^1.5
  rho4 <- kfb(4, gam, lam, tau)/kfb(2, gam, lam, tau)^2
  T <- rho4/8 - 5/24 * rho3^2
  c1 <- 0.5 * log(2) + 0.5 * (p - 1) * log(pi) - 
        0.5 * log( kfb(2, gam, lam, tau) ) - 0.5 * sum( log(lam - tau) ) - 
        tau + 0.25 * sum( gam^2/(lam - tau) )
  ## c1 <- sqrt(2) * pi^(0.5 * (p - 1) ) * kfb(2, gam, lam, tau)^(-0.5) * 
  ## prod(lam - tau)^(-0.5) * exp( -tau + 0.25 * sum( gam^2/(lam - tau) ) )
  c2 <- c1 + log(1 + T)
  c3 <- c1 + T
  ## the next multiplications brings the modification with the negative
  ## values in the lambdas back
  if (mina <= 0) {
    c1 <- c1 + aaa
    c2 <- c2 + aaa
    c3 <- c3 + aaa
  }
  logcon <- c(c1, c2, c3)
  names(logcon) <- c("first order", "second order", "third order")
  logcon
}
