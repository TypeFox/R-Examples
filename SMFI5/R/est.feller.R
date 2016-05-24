est.feller <-
function(  data, method = 'Hessian', days = 360 , significanceLevel = 0.95){
  # Estimation of the parameters of the Feller process
  #
  #        dr = alpha(beta-r)dt + sigma sqrt(r) dW 
  #
  # The time scale is in years and the units are percentages.
  #
  ###########################################################################
  # Input
  #      data: annual bonds yields in percentage;
  #      method: 'Hessian' (default),  'num';
  #      days: number of days per year (default: 360);
  #      significanceLevel: (95\% default).
  #
  # Output
  #       param: parameters (alpha, beta, sigma) of the model;
  #       error: estimation errors for the given confidence level.
  #
  # Example:
    #  out = est.feller(data);
  ###########################################################################
  
  ##
  R   <- data
  h <- 1/days
  
  
  phi0 <- acf(R,1, plot = FALSE)
  alpha0 <- -log(phi0[1]$acf)/h
  beta0  <- mean(R)
  sigma0 <- sd(R)*sqrt(2*alpha0/beta0)
  
  theta0 <- c(log(alpha0), log(beta0), log(sigma0))
  
  
  n <- length( R )
  optim.results <- optim(theta0, function(x) sum(LogLikFeller(x, R, days, n)), hessian=TRUE)
  theta <- optim.results$par
  hessian <- optim.results$hessian
  
  param <- exp(theta)
  
  print(sprintf('\n\nThe estimated coefficients correspond to the annualized spot rate (%%)'))
  
  lmin <- min(eigen(hessian)$values)
  if(lmin<0){
    cat(sprintf('\n The Hessian is not positive definite. Estimation errors are not reliable.\n'))
    cat(sprintf('The numerical method will be used instead.\n\n'))
    method <- 'num'
  }
  
  #  Fisher information matrix 
  
  ##  Fisher information matrix 
  if(method == 'Hessian'){
    FI <- hessian/n
    cat(sprintf('\n Fisher information computed with the numerical Hessian from fminunc (Appendix B.5.1)\n\n'))
  }else{ if(method == 'num'){
    J <- num.jacobian( function(x) LogLikFeller(x, R, days,n),theta,0.01)
    cat(sprintf('\n Fisher information computed with the numerical gradient (Appendix B.5.1)\n\n'))
    FI <- cov(J)
  }}
  
  D <- diag(param)
#  library('corpcor')
  cov_est <- D%*%pseudoinverse(FI)%*%D
  
  ## precision
  criticalValue <- qnorm( 0.5 *(1+significanceLevel) )
  error <- criticalValue * sqrt( diag( cov_est / n ) )
  
  alpha0 <- 0.55
  phi <- exp(-alpha0*h)
  phiest <- exp(-param[1]*h)
  errorphi <- h*phiest*error[1]
  cat(sprintf( '\n alpha <- %.4f /+ %.4f \n', param[1], error[1]))
  cat(sprintf( '\n  beta <- %.4f /+ %.4f \n', param[2], error[2]))
  cat(sprintf( '\n sigma <- %.4f /+ %.4f \n', param[3], error[3]))
  cat(sprintf( '\n   phi <- %.4f phiest <- %.4f /+ %.4f \n', phi, phiest, errorphi))   
  
  return(list(param=param, error=error))
}
