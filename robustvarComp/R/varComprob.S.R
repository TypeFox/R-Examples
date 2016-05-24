#############################################################
# 
#	varComprob.S function
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: June, 24, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

varComprob.S <- function(y, x, V, beta=NULL, gamma=NULL, eta0=NULL, scale=10, control=varComprob.control(), ...) {

# y: matrix. dim(y)=c(p,n) 
# x: array. dim(x)=c(p,n,k)
# V: array. dim(V)=c(p,p,R)
# beta: vector or NULL. length(beta)=k
# gamma: vector or NULL. length(gamma)=R
# eta0: scalar or NULL.
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
  
  initbeta <- beta
  initgamma <- gamma
  initeta0 <- eta0
  initscale <- scale
  
  Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically in Vprod
  v <- qchisq(seq(0.0001,0.9999,length=5000), nrow(Sigma))
  if (control$psi!="rocke")
    s0 <- doSstep(m=v, scale=1, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))
  else
    s0 <- doSsteprocke(m=v, scale=1, bb=control$bb, p=p, arp=control$arp.chi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))

##BEGIN# Iterations
  iter <- 0
  dbeta <- control$rel.tol.beta+1
  dgamma <- control$rel.tol.gamma+1
  dscale <- control$rel.tol.scale+1
  while ((dbeta > control$rel.tol.beta  | dgamma > control$rel.tol.gamma | dscale > control$rel.tol.scale) & iter < control$max.it) {
    iter <- iter+1
    Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically in Vprod
    Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
    Sigmastarinv <- solve(Sigmastar)    
    if (k==0) {
      rr <- y
      beta <- beta1 <- vector(mode="numeric", length=0)
      control$cov <- FALSE
    } else
      rr <- vcrobresid(y=y, x=x, beta=beta)
    RR <- rep(0, n)    
    for (i in 1:n)
      RR[i] <- drop(rr[,i]%*%Sigmastarinv%*%rr[,i])
    if (control$psi!="rocke") {
      scale1 <- doSstep(m=RR, scale=scale, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale/1000, verbose=(control$trace.lev>2))
      W <- vcrobweights(m=RR, scale=scale1, cc=control$tuning.chi, psi=control$psi)
    } else {
      scale1 <- doSsteprocke(m=RR, scale=scale, bb=control$bb, p=p, arp=control$arp.chi, tol=control$rel.tol.scale/1000, verbose=(control$trace.lev>2))
      W <- vcrobweightsrocke(m=RR, scale=scale1, p=p, arp=control$arp.chi)
    }
    Wdot <- scale*W/drop(W%*%RR)
    if (k > 0) {
      XX <- matrix(0, nrow=k, ncol=k)
      XY <- matrix(0, nrow=k, ncol=1)
      for (i in 1:n) {
        XX <- XX + Wdot[i]*t(x[,i,])%*%Sigmastarinv%*%x[,i,] 
        XY <- XY + Wdot[i]*t(x[,i,])%*%Sigmastarinv%*%y[,i]     
      }
      beta1 <- drop(solve(XX)%*%XY)
      dbeta <- max(abs(beta-beta1))
    } else
      dbeta <- 0
    gamma1 <- drop(doGammaClassicSstep(gamma=gamma, resid=rr, scale=scale, V=V, control=control))
    dscale <- max(abs(scale-scale1))
    dgamma <- max(abs(gamma1-gamma))
    if (iter > control$max.it/2) {
      beta <- (beta1+beta)/2
      scale <- (scale1+scale)/2
      gamma <- (gamma1+gamma)/2
    } else {
      beta <- beta1
      scale <- scale1
      gamma <- gamma1
    }    
    if (control$trace.lev>1) {
      cat('Iterations: ', iter, '\n')
      cat('Summary of the weights\n')
      print(summary(c(W)))
      cat('beta: ', beta, '\n')
      cat('gamma: ', gamma, '\n')
      cat('scale: ', scale, '\n')      
      cat('diff max(abs(beta_i - beta_i+1)): ', dbeta, '\n')
      cat('diff max(abs(gamma_i - gamma_i+1)): ', dgamma, '\n')
      cat('diff max(abs(scale_i - scale_i+1)): ', dscale, '\n')
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
  if (control$cov) {
    vcov <- VCOV.ClassicS(beta=beta, gamma=gamma, scale=scale, y=y, x=x, V=V, control=control)
    vcov.beta <- vcov[1:k,1:k]
    vcov.gamma <- vcov[(k+1):(k+R),(k+1):(k+R)]  
  } else {
    vcov.beta <- matrix(NA, k, k)
    vcov.gamma <- matrix(NA, R, R)
  }
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
  result$dotweights <- Wdot
  result$Sigma <- eta0*Vprod(V=V, gamma=gamma)
  result$scale <- result$min <- scale
  result$scale0 <- s0
  result$initial.values <- list()  
  result$initial.values$beta <- initbeta
  result$initial.values$gamma <- initgamma
  result$initial.values$eta0 <- initeta0
  result$initial.values$scale <- initscale  
  result$iterations <- iter
  result$control <- control
  result$control$method <- "S"  
  class(result) <- 'varComprob.S'
  return(result)
}
