#Source: artSim.R
artSim <- function(obj=NULL, n=100, d=0, lambda=0, phi=numeric(0), 
                   theta=numeric(0), mean=0, sigma2=1) {
    if (identical(class(obj),"ARTFIMA")) {
      d <- obj$dHat
      lambda <- obj$lambdaHat
      mean <- obj$constant
      sigma2 <- obj$sigmaSq
      phi <- obj$phi
      theta <- obj$theta
    }
  r <- tacvfARTFIMA(d=d, lambda=lambda, phi=phi, theta=theta, maxlag=n-1)
  mean+sigma2*DLSimulate(n, r)
}