LogLikVasicek <-
function( theta, R, tau, days, n ){
  #  This function computes the -log-likelihoods for the Vasicek model
  #
  #           dr = alpha(beta-r)dt + sigma  dW,
  #
  # with market price of risk q(r) = q1+q2 r.
  #
  # Input
  #   theta: parameters of the annualized spot rate in percent;
  #   R: annual returns in percent of the bonds (n x 1);
  #   tau: maturities (n x 1) in years;
  #   scalingFact: measurement scale ( 100 x h, h is 1 if annual, 360 or 365
                                       #                if daily).  
  #
  # Output
  #    -log-likelihoods (LL).
  
  
  scalingFact <- 100
  likelihoodFactor <- 0.5 * log( 2 * pi )
  h <- 1/days
  
  LL <- matrix(1,n-1,1)
  
  alpha <- exp(theta[1])
  beta  <- exp(theta[2])
  sigma <- exp(theta[3])
  q1    <- theta[4]
  q2    <- theta[5]
  
  param <- c(alpha, beta, sigma, q1, q2)
  
  
  params <- get.vasicek.param(param, tau, scalingFact )
  A <- params$A
  B <- params$B
  if((params$a<=0) || (params$b<=0)){
    LL <- 1.0e+20*LL
    return(LL)
    
  }
  
  r <- ( tau * R + scalingFact * A ) / B
  if(sum( r <= 0 ) >0){
    LL <- 1.0e+20*LL
    return(LL)
  }
  
  # parameters for the returns\
  phi     <- exp( - alpha *h)
  gamma   <- sigma * sqrt( ( 1 - phi^2 ) / ( 2 * alpha ) )
  
  eps <- ( r[-1] - beta - phi * ( head(r,length(r)-1) - beta ) ) / gamma
  
  # log-likelihood function to be minimized
  LL = likelihoodFactor + 0.5 * eps^2 + log( gamma ) + log(  B[-1] / tau[-1]  )
  
  if(any(is.infinite( LL )) || any(is.nan( LL )) || any(is.complex( LL ))){
    LL <- 1.0e+20*matrix(1,n-1,1)
    return(LL)
  
  }
  return(LL)
}
