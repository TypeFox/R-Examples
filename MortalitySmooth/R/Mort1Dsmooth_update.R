Mort1Dsmooth_update <-
function(x, y, offset, wei, psi2, B,
                                lambdaP, a){
  ## Input:
  ## x: abcissae of data
  ## y: count response
  ## offset: an a priori known component
  ## wei: weights
  ## psi2: overdispersion parameter
  ## B: B-splines basis
  ## P: penalty term/matrix
  ## a: coefficients
  
  ## Output:
  ## a: updated coefficients
  
  ## linear predictor
  eta <- B%*%a
  ## expected values
  mu <- exp(eta + offset)
  ## weights
  w <- c(wei*mu)
  ## working response
  z <- c(wei*((y - mu)/mu + eta))
  ## regression
  BtWB <- t(B) %*% (w * B)
  BtWz <- t(B) %*% (w * z)
  a <- solve(BtWB + psi2*lambdaP, BtWz)
  ## coefficients
  a <- matrix(a, ncol = 1)
  a
}
