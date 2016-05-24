LogLikOU <-
function( theta, R, days, n ){
  #  This function computes the -log-likelihoods for the Orstein-Uhlenbeck model
  #
  #dr = alpha(beta-r)dt + sigma  dW.
  #
  # Input
  #      theta: log of parameters of the annualized spot rate; 
  #           R: annual returns in percent of the bonds (n x 1);
  #        days: number of days in a year. 
  #
  # Output
  #        LL: -log-likelihoods.
  
  
  h <- 1/days
  likelihoodFactor <- 0.5 * log( 2 * pi )
  
  LL <- matrix(1,n-1,1)
  
  alpha <- exp(theta[1])
  beta  <- theta[2]
  sigma <- exp(theta[3])
  
  phi     <- exp( - alpha *h)
  gamma   <- sigma * sqrt( ( 1 - phi^2 ) / ( 2 * alpha ) )
  
  eps <- ( R[-1] - beta - phi * ( head(R,length(R)-1) - beta ) ) / gamma
  
  # log-likelihood function to be minimized
  LL = likelihoodFactor + 0.5 * eps^2 + log( gamma )
  
  if(any(is.infinite( LL )) || any(is.nan( LL )) || any(is.complex( LL ))){
    LL <- 1.0e+20*matrix(1,n-1,1)
    return(LL)
  }
  
  return(LL)
  
}
