LogLikCIR <-
function( theta, R, tau, days, n ){
  #  This function computes the -Log-likelihoods for the CIR model
  #
  #            dr <- alpha(beta-r)dt + sigma sqrt(r) dW   
  # with market price of risk q(r) <- q1/sqrt(r) +q2 sqrt(r)
  #
  # Input
  #       theta: parameters of the annualized spot rate 
  #           R: annual returns in percent of the bonds (n x 1)
  #         tau: maturities (n x 1) in days
  #        days: number of days in a year.
  #
  # Output
  #        LL: -log-likelihoods.
  

  
  scalingFact <- 100
  h <- 1/days
  
  LL <- matrix(1,n-1,1)
  
  alpha <- exp(theta[1])
  beta  <- exp(theta[2])
  sigma <- exp(theta[3])
  q1    <- theta[4]
  q2    <- theta[5]
  
  param <- c(alpha, beta, sigma, q1, q2)
  
  
  params <- get.cir.param(param, tau, scalingFact )
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
  
  # parameters for the returns
  phi     <- exp( - alpha *h)
  nu      <- 4 * alpha * beta / sigma^2
  omega <- beta * ( 1 - phi ) / nu  
  
  r       <- r / omega
  D       <- head(r, n=-1) * phi
  
  # log-likelihood function
  #LL <- log( omega ) +  log( abs( B(2:}) ./ tau(2:}) ) ) -  log( max(  PDFNChi2( r(2:}), nu, D ), eps ) ) 
  # Better take the Matlab function ncx2pdf which is more stable numerically.
  LL <- log( omega ) +  log( abs( B[-1] / tau[-1] ))  -  log( pmax(  dchisq( r[-1], nu, D ), .Machine$double.ep ) ) 
  
  if(any(is.infinite( LL )) || any(is.nan( LL )) || any(is.complex( LL ))){
    LL <- 1.0e+20*matrix(1,n-1,1)
    return(LL)
  }
  
  return(LL)
}
