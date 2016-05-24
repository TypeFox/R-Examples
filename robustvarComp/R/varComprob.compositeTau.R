#############################################################
# 
#	varComprob.compositeTau function
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: June, 28, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

varComprob.compositeTau <- function(y, x, V, beta=NULL, gamma=NULL, eta0=NULL, scales=rep(10, p*(p-1)/2), control=varComprob.control(), ...) {
  G  <- 3.25 # for rho=optimal, see for instance optimChi.f
# y: matrix. dim(y)=c(p,n) 
# x: array. dim(x)=c(p,n,k)
# V: array. dim(V)=c(p,p,R)
# beta: vector or NULL. length(beta)=k
# gamma: vector or NULL. length(gamma)=R
# eta0: scalar or NULL.  
# Sigma: matrix or NULL. dim(Sigma)=c(p,p)  
# scales: vector or NULL. length(scales)=p*(p-1)/2

  if (control$psi=="rocke")
    stop("Rocke rho function is not yet available for composite methods")
  
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
  initscales <- scales
  
  v <- qchisq(seq(0.0001,0.9999,length=5000), 2)
  s0 <- doSstep(m=v, scale=1, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))

##BEGIN# Iterations
  iter <- 0
  dbeta <- control$rel.tol.beta+1
  dgamma <- control$rel.tol.gamma+1
  dscale <- control$rel.tol.scale+1
  while ((dbeta > control$rel.tol.beta  | dgamma > control$rel.tol.gamma | dscale > control$rel.tol.scale) & iter < control$max.it) {
    iter <- iter+1
    Sigma <- Vprod(V=V, gamma=gamma) ## V0 is added automatically in Vprod
    if (k==0) {
      rr <- y
      beta <- beta1 <- vector(mode="numeric", length=0)
      control$cov <- FALSE
    } else
      rr <- vcrobresid(y=y, x=x, beta=beta)
    RR <- rssr(resid=rr, Sigma=Sigma)
    scales1 <- doSsteppw(RR=RR, scale=scales, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale/1000, verbose=(control$trace.lev>2))
    if (control$psi=="optimal")
      scales2 <- scales1*control$tuning.chi^2
    else
      scales2 <- scales1
    tau <- doTausteppw(RR=RR, scale=scales2, cc=control$tuning.psi, psi=control$psi)
    W <- vcrobweightspw(RR=RR, scale=scales1, cc=control$tuning.chi, psi=control$psi)
    W2 <- vcrobweightspw(RR=RR, scale=scales2, cc=control$tuning.psi, psi=control$psi)    
    Wdot <- vcrobweightsdotpw(W=W, RR=RR, scale=scales2)
    rhospw <- rhostarpw(RR=RR, scale=scales2, cc=control$tuning.psi, psi=control$psi)
    WW <- vcrobweights2pw(RR=RR, scale=scales2, cc=control$tuning.psi, psi=control$psi)
    ### These modifications are needed because the rho is normalized in [0,1]
    if (control$psi=="bisquare")
      Wdot <- Wdot*(1/3*control$tuning.psi^2*rhospw - WW)
    else if (control$psi=="optimal")
      Wdot <- Wdot*(2*G*control$tuning.psi^2*rhospw - WW)
    if (k > 0) {
      XX <- xssx(x=x, Sigma=Sigma)
      XY <- xssy(x=x, y=y, Sigma=Sigma)
      beta1 <- drop(doBetastep(Wdot=(Wdot+W2), XX=XX, XY=XY))
      dbeta <- max(abs(beta-beta1))
    } else
      dbeta <- 0
    gamma1 <- drop(doGammaCompositeTaustep(gamma=gamma, resid=rr, scales=scales1, V=V, Tmax=tau+10, control=control))
    dscale <- max(abs(scales-scales1))
    dgamma <- max(abs(gamma1-gamma))
    if (iter > control$max.it/2) {
      beta <- (beta1+beta)/2
      scales <- (scales1+scales)/2
      gamma <- (gamma1+gamma)/2
    } else {
      beta <- beta1
      scales <- scales1
      gamma <- gamma1
    }    
    if (control$trace.lev>1) {
      cat('Iterations: ', iter, '\n')
      cat('Summary of the weights\n')
      print(summary(c(W)))
      cat('beta: ', beta, '\n')
      cat('gamma: ', gamma, '\n')
      cat('sum(scales): ', sum(scales), '\n')
      cat('Tau: ', tau, '\n')
      cat('diff max(abs(beta_i - beta_i+1)): ', dbeta, '\n')
      cat('diff max(abs(gamma_i - gamma_i+1)): ', dgamma, '\n')
      cat('diff max(abs(scale_i - scale_i+1)): ', dscale, '\n')
    }
  }
##END# Iterations

##BEGIN# Eta0
  Sigma <- Vprod(V=V, gamma=gamma)
  RSR <- rsr(resid=rr, Sigma=Sigma)  
  eta0 <- doSstep(m=RSR/s0, scale=1, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale, verbose=(control$trace.lev>2))
##END# Eta0

##BEGIN# VCOV
  if (control$cov) {
    vcov <- VCOV.CompositeTau(beta=beta, gamma=gamma, scales=scales, y=y, x=x, V=V, control=control)
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
  result$weights <- result$weights.1 <- W
  result$weights.2 <- W2  
  result$dotweights <- Wdot
  result$Sigma <- eta0*Vprod(V=V, gamma=gamma)
  result$scales <- scales
  result$tau <- result$min <- tau
  result$scale0 <- s0
  result$initial.values <- list()  
  result$initial.values$beta <- initbeta
  result$initial.values$gamma <- initgamma
  result$initial.values$eta0 <- initeta0
  result$initial.values$scales <- initscales  
  result$iterations <- iter 
  result$control <- control
  result$control$method <- "compositeTau"  
  class(result) <- 'varComprob.compositeTau'
  return(result)
}
