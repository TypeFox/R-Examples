get.cir.param <-
function(param, tau, scalingFact = 1 ){
  # This function computes the terms A and B for the price of a zero-coupon
  # bond under the CIR model.
  
  sqrtScalingFact <- sqrt( scalingFact )
  
  alpha <- param[1]
  beta  <- param[2] / scalingFact
  sigma <- param[3] / sqrtScalingFact
  q1    <- param[4] / sqrtScalingFact
  q2    <- param[5] * sqrtScalingFact
  
  a <- alpha + q2 * sigma
  
  b <- ( alpha * beta - q1 * sigma ) / a
  
  gam <- sqrt( a^2 + 2 * sigma^2 )
  
  d <- 1 - exp( - gam * tau )
  
  B <- d / ( 0.5 * ( gam + a ) * d + gam * ( 1 - d ) )
  
  A <- - ( 2 * a * b / sigma^2 ) * ( 0.5 * ( gam - a ) * tau + log( 0.5 * ( 1 + a / gam ) * d + ( 1 - d ) ) ) 
  
  return(list(A=A, B=B, a=a, b=b))
  
}
