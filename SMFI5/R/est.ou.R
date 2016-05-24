est.ou <-
function(  data, method = 'Hessian', days = 360 , significanceLevel = 0.95){
  #  Estimation of the parameters of the OU process
  #
  #           dr = alpha(beta-r)dt + sigma  dW,
  #
  #
  #
  ###########################################################################
  # Input
  #      data: observations;
  #      method: 'Hessian' (default), 'num';
  #      days: number of pdays per year (default: 360);
  #      significanceLevel: (95# default).
  #
  # Output
  #       paramNum: Numerical estimation of the parameters alpha, beta, sigma, phi;
  #       errorNum: Estimation errors for a 95# confidence level of paramNum;
  #       paramExp: Explicit estimation of the parameters alpha, beta, sigma, phi;
  #       errorExp: Estimation errors for a 95# confidence level of paramExp.
  #
  # Example:
    #  out = est.OU(data);
  ###########################################################################
  
  ##
  R   <- data
  h <- 1/days
  
  theta0 <- c(log(1), mean(R), log(1))
  
  
  n <- length( R )
  optim.results <- optim(theta0, function(x) sum(LogLikOU(x, R, days, n)), hessian=TRUE)
  theta <- optim.results$par
  hessian <- optim.results$hessian
  
  alpha <- exp(theta[1])
  beta  <- theta[2]
  sigma <- exp(theta[3])
  phi <- exp(-alpha * h )
  
  paramNum <- c(alpha,beta,sigma,phi)
  
  
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
    J <- num.jacobian( function(x) LogLikOU(x, R, days,n),theta,0.01)
    cat(sprintf('\n Fisher information computed with the numerical gradient (Appendix B.5.1)\n\n'))
    FI <- cov(J)
  }}
  
  D <- diag(c(alpha, 1, sigma))
 # library('corpcor')
  cov_est <- D%*%pseudoinverse(FI)%*%D
  
  ## precision
  criticalValue <- qnorm( 0.5 *(1+significanceLevel) )
  errorNum0 <- criticalValue * sqrt( diag( cov_est / n ) )
  errorNum4 <- phi * h * errorNum0[1]
  errorNum <- c(errorNum0, errorNum4)
  
  cat(sprintf( '\n alpha <- %.4f /+ %.4f \n', paramNum[1], errorNum[1]))
  cat(sprintf( '\n  beta <- %.4f /+ %.4f \n', paramNum[2], errorNum[2]))
  cat(sprintf( '\n sigma <- %.4f /+ %.4f \n', paramNum[3], errorNum[3]))
  cat(sprintf( '\n   phi <- %.4f /+ %.4f \n', paramNum[3], errorNum[4]))   
  
  
  ## exp MLE estimation
  
  z1 <- R[-1]
  z2 <- head(R,length(R)-1)
  m1 <- mean(z1)
  m2 <- mean(z2)
  s1 <- sd(z1) 
  s2 <- sd(z2)
  
  h <- 1/days
  rho <- cor(z1,z2)
  
  phiest <- rho*s1/s2
  
  gamma <- sd( z1-phiest*z2)
  g <- 1-phiest^2
  
  alphaexp <- -log(phiest)/h
  betaexp <- (m1-phiest*m2)/(1-phiest)
  sigmaexp <- gamma*sqrt(2*alphaexp/g)
  zeta <- (sigmaexp^2/gamma)*( 2*alphaexp*phiest^2*h-g  )/(phiest*g^2*h )
  
  
  error_exp <- mat.or.vec(5,1)
  error_exp[4] <- criticalValue*sqrt(g/n)
  error_exp[1] <- error_exp[4]/(h*phiest)
  error_exp[2] <- criticalValue*gamma/(1-phiest)/sqrt(n)
  error_exp[3] <- criticalValue*sqrt((zeta^2*g + 0.5*sigmaexp^2)/n)
  error_exp[5] <- criticalValue*sqrt(0.5*gamma^2/n)
  
  
  
  Vexp <- diag(c(g,gamma^2/(1-phiest)^2,0.5*gamma^2))
  
  Jexp <- matrix(c( -1/(phiest*h), 0, 0, 0, 1, 0, zeta, 0, sigmaexp/gamma,  1, 0, 0, 0, 0, 1), ncol=3, nrow=5)
  
  covExp <- Jexp%*%Vexp%*%t(Jexp)
     
  paramExp <- c(alphaexp,betaexp,sigmaexp,phiest,gamma)
  errorExp <- criticalValue*sqrt(diag(covExp)/n)
  
  cat(sprintf('\n\n Estimation with the explicit method\n\n'))
  
  
  cat(sprintf( '\n alpha = %.4f /+ %.4f \n', alphaexp, error_exp[1]))
  cat(sprintf( '\n  beta = %.4f /+ %.4f \n', betaexp, error_exp[2]))
  cat(sprintf( '\n sigma = %.4f /+ %.4f \n', sigmaexp, error_exp[3]))
  cat(sprintf( '\n phi   = %.4f /+ %.4f \n', phiest, error_exp[4]))
  cat(sprintf( '\n gamma = %.4f /+ %.4f \n', gamma, error_exp[5]))
  
  
  
  
  
  return(list(param.num=paramNum, error.num=errorNum, param.exp=paramExp, error.exp=errorExp))
}
