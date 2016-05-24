LogLikFeller <-
function( theta, R, days, n ){
  #  This function computes the -Log-likelihoods for the CIR model
  #
  # dr = alpha(beta-r)dt + sigma sqrt(r) dW   
  #
  # Input
  #      theta: log of parameters of the annualized spot rate; 
  #           R: annual returns in percent of the bonds (n x 1);
  #        days: number of days in a year. 
  #
  # Output
  #        LL: -log-likelihoods.
  
  
  h <- 1/days
  
  LL <- matrix(1,n-1,1)
  
  alpha <- exp(theta[1])
  beta  <- exp(theta[2])
  sigma <- exp(theta[3])

  phi     <- exp( - alpha *h)
  nu      <- 4 * alpha * beta / sigma^2
  omega   <- beta * ( 1 - phi )/ nu
  z       <- R / omega
  D       <- head(z,length(z)-1) * phi # non-centrality parameters
  
  # log-likelihood function to be minimized
  LL <- log( omega ) - log( pmax(  dchisq( z[-1], nu, D ), .Machine$double.ep ) )
  
  if(any(is.infinite( LL )) || any(is.nan( LL )) || any(is.complex( LL ))){
    LL <- 1.0e+20*matrix(1,n-1,1)
    return(LL)
  }
  
  return(LL)
  
}
