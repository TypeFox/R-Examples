est.cir <-
function(  data, method = 'Hessian', days = 360 , significanceLevel = 0.95){
  # Estimation of the parameters of the CIR model
  #
  #        dr <- alpha(beta-r)dt + sigma sqrt(r) dW 
  #
  # with market price of risk q(r) <- q1/sqrt(r) +q2 sqrt(r). The time scale 
  # is in years and the units are percentages.
  #
  ###########################################################################
  # Input
  #      data: [R,tau] (n x 2), with R: annual bonds yields in percentage,
  #             and tau: maturities in years
  #      method: 'Hessian' (default), 'num'
  #      days: number of days per year (default: 360)
  #      significanceLevel: (95# default).
  #
  # Output
  #       param: parameters (alpha, beta, sigma, q1,q2) of the model
  #       error: estimation errors for the given confidence level
  #       rimp: implied spot rate.
  #
  # Example:
  #  out <- est.cir(data.cir)
  ###########################################################################
  
  
  ##
  R   <- data[,1]
  tau <- data[,2]
  h <- 1/days
  scalingFact <- 100
  
  ## starting values
  
  # Rough estimation in the case of a Feller process
  phi0 <- acf(R,1, plot = FALSE)
  alpha0 <- -log(phi0[1]$acf)/h
  beta0  <- mean(R)
  sigma0 <- sd(R)*sqrt(2*alpha0/beta0)
  nu0 <- 0.25*sigma0^2/(alpha0*beta0)
  param0 <-  c(alpha0, beta0, sigma0, 0 ,0)
  theta0 <- c(log(alpha0), log(beta0), log(sigma0), 0,  0)
  #theta0 <- [log(2) log(mean(R)) log(1) 0 0]
  #param0 <- [0.52.550.3650.30] # true parameters
  
  
  # pricing coefficients (computed by changing the scale from percentage
  # to fraction)
  
  params <- get.cir.param( param0, tau, scalingFact )
  A <- params$A
  B <- params$B
  # get short rate
  r <- ( tau * R + scalingFact * A )/ B    
  
  plot(r,,'l')
  title('Annual implied rates (%) for the starting parameters')
  
  ANSWER <- readline("Are you a satisfied with the graph? ")
  
  if (substr(ANSWER, 1, 1) == "n"){
    cat("\nIt the rates are too big, decrease q1.\n")
    cat("\nIf the rates are too small, increase q1.\n")
    cat(sprintf('\n\n'))}
  
  else
  {
#   alpha   <- NaN
#   beta    <- NaN
#   sigma   <- NaN
#   q1      <- NaN
#   q2      <- NaN
#   error   <- NaN
#   r       <- NaN
  
  ## Options for fminunc
  # options <- optimset( 'largescale',   'on', ...
  #                     'MaxFunEvals',  50000, ...
  #                     'MaxIter',      5000, ...
  #                     'Display',      'iter', ...
  #                     'TolFun'      , 1e-6, ...
  #                     'TolCon'      , 1e-6, ...
  #                     'TolX'        , 1e-6, ...                        
  #                     'diagnostics',  'on')
  
  ## Maximisation of the log-likelihood
  n <- length( R )
  optim.results <- optim(theta0, function(x) sum(LogLikCIR(x, R, tau, days, n)), hessian=TRUE)
  theta <- optim.results$par
  hessian <- optim.results$hessian
  alpha <- exp(theta[1])
  beta  <- exp(theta[2])
  sigma <- exp(theta[3])
  q1    <- theta[4]
  q2    <- theta[5]
  
  param <- c(alpha,beta,sigma,q1,q2)
  cir.param.est <- get.cir.param(param, tau, scalingFact )
  A<- cir.param.est$A
  B <- cir.param.est$B
  
  # Implied spot rates
  r <- ( tau * R + scalingFact * A ) / B
  
  # real parameters
  paramTrue <- c(0.5,2.55,0.365,0.3,0)
  
  
  cir.param.true <- get.cir.param( paramTrue, tau, scalingFact )
  
  A <- cir.param.true$A
  B <- cir.param.true$B
  # Implied spot rates
  r0 <- ( tau * R + scalingFact * A ) / B
  
  
  cat(sprintf('\n\nThe estimated coefficients correspond to the annualized spot rate (##)'))
  
  
 # library('ggplot2')
 # library('reshape')
  tmp <- data.frame(x=(1:length(r)), r=r, r0=r0)
  names(tmp) <- c('x',c('Implied rates', 'Real rates'))
  framed.data = melt(tmp, id='x')
  
  # To remove the notes by R CMD check!
  value <- NULL 
  variable <- NULL
  x <- NULL
  
  values.graph <- ggplot(framed.data,
                         aes(x=x, y=value, group=variable, colour=variable))
  
  values.graph <- values.graph + geom_line() + 
    ggtitle('Implied annual spot rates (#)')
  
  print(values.graph)
  
  n <- length( R )
  
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
    J <- num.jacobian( function(x) LogLikCIR(x, R, tau, days,n),theta,0.01)
    cat(sprintf('\n Fisher information computed with the numerical gradient (Appendix B.5.1)\n\n'))
    FI <- cov(J)
  }}
  
  D <- diag(c(alpha,beta,sigma,1,1))
 # library('corpcor')
  cov_est <- D%*%pseudoinverse(FI)%*%D
  
  ## precision
  criticalValue <- qnorm( 0.5 *(1+significanceLevel) )
  error <- criticalValue * sqrt( diag( cov_est / n ) )
  
  phi <- exp(-paramTrue[1]*h)
  phiest <- exp(-param[1]*h)
  errorphi <- h*phiest*error[1]
  cat(sprintf( '\n alpha = %.4f /+ %.4f \n', param[1], error[1]))
  cat(sprintf( '\n  beta = %.4f /+ %.4f \n', param[2], error[2]))
  cat(sprintf( '\n sigma = %.4f /+ %.4f \n', param[3], error[3]))
  cat(sprintf( '\n    q1 = %.4f /+ %.4f \n', param[4], error[4]))
  cat(sprintf( '\n    q2 = %.4f /+ %.4f \n', param[5], error[5]))
  cat(sprintf( '\n   phi = %.4f,  phiest = %.4f /+ %.4f \n', phi, phiest, errorphi))   
  
  return(list(param=param, error=error, r=r))
  }
}
