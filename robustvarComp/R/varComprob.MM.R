#############################################################
# 
#	varComprob.MM function
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: June 24, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

varComprob.MM <- function(y, x, V, S=NULL, beta=NULL, gamma=NULL, scale=NULL, control=varComprob.control()) {

# y: matrix. dim(y)=c(p,n) 
# x: array. dim(x)=c(p,n,k)
# V: array. dim(V)=c(p,p,R)  
# beta: vector or NULL. length(beta)=k
# gamma: vector or NULL. length(gamma)=R
# scale: scalar or NULL.
  
## SET STORAGE MODE OF y, x and V
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(V) <- "double"

  xdim <- dim(x)
  p <- xdim[1]
  n <- xdim[2]
  k <- xdim[3]
  Vdim <- dim(V)
  R <- Vdim[3]
  JL <- p*(p-1)/2

  
## About initial values  
  if (is.null(beta))
    beta <- S$beta
  if (is.null(gamma))
    gamma <- S$gamma
  if (is.null(scale))
    scale <- S$scale
  
  if(is.null(beta) | is.null(gamma) | is.null(scale))
    stop('Initial values for beta, gamma and value for the scale must be supplied')

  v <- qchisq(seq(0.0001,0.9999,length=5000), p)
  if (control$psi!="rocke")
    s0 <- doSstep(m=v, scale=1, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))
  else
    s0 <- doSsteprocke(m=v, scale=1, bb=control$bb, p=p, arp=control$arp.chi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))
  
##BEGIN# Iterations
  iter <- 0
  dbeta <- control$rel.tol.beta+1
  dgamma <- control$rel.tol.gamma+1  
  while ((max(dbeta) > control$rel.tol.beta | dgamma > control$rel.tol.gamma)  & iter < control$max.it) {
    iter <- iter+1
    ## SIGMA
    Sigma <- Vprod(V=V, gamma=gamma)
    Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
    Sigmastarinv <- solve(Sigmastar)        
    ## RESIDUALS
    if (k==0) {
      rr <- y
      beta <- beta1 <- vector(mode="numeric", length=0)
      control$cov <- FALSE
    } else
      rr <- vcrobresid(y=y, x=x, beta=beta)
    ## SQUARED MAHALANOBIS DISTANCES
    RR <- rep(0, n)
    for (i in 1:n)
      RR[i] <- drop(rr[,i]%*%Sigmastarinv%*%rr[,i])
    ## WEIGHTS
    if (control$psi!="rocke")
      W <- vcrobweights(m=RR, scale=scale, cc=control$tuning.psi, psi=control$psi)
    else
      W <- vcrobweightsrocke(m=RR, scale=scale, p=p, arp=control$arp.psi)
    ## SQUARED MAHALANOBIS DISTANCES
    if (k > 0) {
      XX <- matrix(0, nrow=k, ncol=k)
      XY <- matrix(0, nrow=k, ncol=1)    
      for (i in 1:n) {      
        XX <- XX + W[i]*t(x[,i,])%*%Sigmastarinv%*%x[,i,] 
        XY <- XY + W[i]*t(x[,i,])%*%Sigmastarinv%*%y[,i]     
      }
      ## BETAS
      beta1 <- drop(solve(XX)%*%XY)
      dbeta <- max(abs(beta-beta1))
    } else
      dbeta <- 0
    ## GAMMAS
    Mmax <- doGammaClassicMMGoal(x=gamma, resid=rr, scale=scale, V=V, Mmax=NA, controllo=control)+10
    if (is.na(Mmax))
      stop("The Sigma matrix is singular and we do not know how to fix it")
    gamma1 <- drop(doGammaClassicMMstep(gamma=gamma, resid=rr, scale=scale, V=V, Mmax=Mmax, control=control))
    dgamma <- max(abs(gamma1-gamma))
    if (iter > control$max.it/2) {
      beta <- (beta1+beta)/2
      gamma <- (gamma1+gamma)/2
    } else {
      beta <- beta1
      gamma <- gamma1
    }
    if (control$trace.lev>1) {
      cat('Iterations: ', iter, '\n')
      cat('beta: ', beta, '\n')
      cat('gamma: ', gamma, '\n')
      M <- doGammaClassicMMGoal(x=gamma, resid=rr, scale=scale, V=V, Mmax=Mmax, controllo=control)
      cat('M: ', M, '\n')
      cat('diff max(abs(beta_i - beta_i+1)): ', dbeta, '\n')
      cat('diff max(abs(gamma_i - gamma_i+1)): ', dgamma, '\n')
    }
  }
##END# Iterations
  
##BEGIN# Eta0
  Sigma <- Vprod(V=V, gamma=gamma)
  Sigmainv <- solve(Sigma)
  RSR <- rep(0, ncol(rr))  
  for (i in 1:ncol(rr))
    RSR[i] <- drop(rr[,i]%*%Sigmainv%*%rr[,i])
  if (control$psi!="rocke")
    eta0 <- doSstep(m=RSR/s0, scale=1, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))
  else
    eta0 <- doSsteprocke(m=RSR/s0, scale=1, bb=control$bb, p=p, arp=control$arp.chi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))
##END# Eta0
  
##BEGIN# VCOV
  vcov <- VCOV.ClassicMM(beta=beta, gamma=gamma, scale=scale, y=y, x=x, V=V, control=control)
  vcov.beta <- vcov[1:k,1:k]
  vcov.gamma <- vcov[(k+1):(k+R),(k+1):(k+R)]  
##END# VCOV

  result <- list()
  result$call <- match.call()
  result$beta <- drop(beta)
  result$vcov.beta <- vcov.beta
  result$eta <- drop(gamma*eta0)  
  result$vcov.eta <- vcov.gamma*eta0^2
  result$gamma <- drop(gamma)
  result$vcov.gamma <- vcov.gamma
  result$eta0 <- eta0    
  result$resid <- rr
  result$weights <- W
  result$scales <- scale
  result$scale0 <- s0
  result$min <- doGammaClassicMMGoal(x=gamma, resid=rr, scale=scale, V=V, Mmax=NA, controllo=control)
  result$iterations <- iter   
  result$control <- control
  result$control$method <- "MM"
  class(result) <- 'varComprob.MM'
  return(result)
}
