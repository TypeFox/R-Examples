context("LagEstimator")

test_that("lagEstimator works as expected",{
  
  set.seed(9247)
  
  Bn = 10
  Y <- rnorm(128)
  
  freq = 2*pi*(0:8)/16
  
  levels.1 <- c(0.1, 0.5, 0.9)
  
  weight <- lagKernelWeight(W = WParzen,  bw = Bn, K = length(Y))
  lagOp <- clippedCov(Y, levels.1 = levels.1)
  lagEst <- lagEstimator(lagOp,weight = weight)
  V <- getValues(lagEst, frequencies = freq)[,,,1]
  
  K <- length(levels.1)
  n <- length(Y)
  
  #Calculate lag-window estimator by hand:
  
  # weight function
  K_Parzen <- function (u) {
    if (abs(u) <= 1) {
      if(abs(u)<=0.5){return(1-6*u^2+6*abs(u)^3)}    
      else{return(2*(1-abs(u))^3)}
      
    } else {
      return(0)
    }
  }
  
  Ker <- Vectorize(K_Parzen)
  
  # We trust that clippedCov works as expected!
  VcC <- getValues(lagOp)[,,,1]
    
  # Define and fill a result vector:
  
  res <- array(dim = c(length(freq), K, K))
  
  ii <- complex(real = 0, imaginary = 1)
  
  for (k1 in 1:K) {
    for (k2 in 1:K) {
      for (f in seq(freq)) {
        # sum for k >= 0
        aux1 <- sum(VcC[1:n, k1, k2] * Ker( (0:(n-1)) / Bn ) * exp(- freq[f] * ii * (0:(n-1))) )
        aux2 <- sum(VcC[2:n, k2, k1] * Ker( -(1:(n-1)) / Bn ) * exp(- freq[f] * ii * -(1:(n-1))) )
        res[f, k1, k2] <- (aux1 + aux2) / (2*pi)
      }
    }
  }
  
  expect_equal(dim(V), dim(res))
  expect_equal(V, res)
}
)