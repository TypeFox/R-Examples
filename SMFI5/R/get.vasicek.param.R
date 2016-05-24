get.vasicek.param <-
function(param, tau, scalingFact = 1){
  # This function computes the terms A and B for the price of a zero-coupon
  # bond under the Vasicek model.
  
  alpha <- param[1]
  beta  <- param[2] / scalingFact
  sigma <- param[3] / scalingFact
  q1    <- param[4]
  q2    <- param[5] * scalingFact
  
  a <- alpha + q2 * sigma
  b <- ( alpha * beta - q1 * sigma) / a
  
  eps <- .001
  if(a > eps){
  B <- ( 1 - exp( - a * tau ) ) / a
  }else{
    B <- tau
  }
  
  A <- - ( b - sigma^2 / ( 2 * a^2 ) ) * ( tau - B ) - sigma^2 * B^2 / ( 4 * a )
  
  return(list(A=A, B=B, a=a, b=b))
  }
