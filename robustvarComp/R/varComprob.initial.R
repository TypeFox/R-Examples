#############################################################
# 
#	varComprob.initial function
#	Author: Claudio Agostinelli and Victor J. Yohai
#	E-mail: claudio@unive.it
#	Date: July, 01, 2014
#	Version: 0.1
#
#	Copyright (C) 2014 Claudio Agostinelli
#                      and Victor J. Yohai
#
#############################################################

varComprob.initial <- function(y, x, V, beta=NULL, gamma=NULL, eta0=NULL, Sigma=NULL, scales=NULL, scale=NULL, control=varComprob.control(), ...) {

# y: matrix. dim(y)=c(p,n) 
# x: array. dim(x)=c(p,n,k)
# V: array. dim(V)=c(p,p,R)
# beta: vector or NULL. length(beta)=k
# gamma: vector or NULL. length(gamma)=R
# eta0: scalar or NULL.  

# Sigma: matrix or NULL. dim(Sigma)=c(p,p)  
# scales: vector or NULL. length(scales)=p*(p-1)/2

## SET STORAGE MODE OF y, x and V
  storage.mode(y) <- "double"
  storage.mode(x) <- "double"
  storage.mode(V) <- "double"

##BEGIN# Initial values
  xdim <- dim(x)
  p <- xdim[1]
  n <- xdim[2]
  k <- xdim[3]
  Vdim <- dim(V)
  R <- Vdim[3]
  JL <- p*(p-1)/2
  if (is.null(beta) & k > 0) {
    xx <- matrix(0, nrow=n*p, ncol=k)
    for (i in 1:k) {
      xx[,i] <- as.vector(x[,,i])
    }
    if (control$beta.univ) {
      beta <- rep(0, k)
      for (i in 1:k) {
        if (control$fixed.init=="lmrob.S")
          beta[i] <- lmrob.S(x=xx[,i], y=as.vector(y), control=lmrob.control(), trace.lev=control$trace.lev)$coef
        else
          beta[i] <- lmRob(I(as.vector(y))~I(xx[,i]) - 1, control=lmRob.control(efficiency=0.85, initial.alg="Fast", final.alg="MM"))$coef
      }
    } else {
      if (control$fixed.init=="lmrob.S")
        beta <- lmrob.S(x=xx, y=as.vector(y), control=lmrob.control(), trace.lev=control$trace.lev)$coef
      else
        beta <- lmRob(I(as.vector(y))~I(xx) - 1, control=lmRob.control(efficiency=0.85, initial.alg="Fast", final.alg="MM"))$coef
    }
  }
  if (k==0)
    rr <- y
  else
    rr <- vcrobresid(y=y, x=x, beta=beta) ## residuals 
  if (is.null(gamma)) {
    if (is.null(Sigma)) {
      if (is.function(control$cov.init))
        Sigma <- control$cov.init(t(rr), ...)
      else if (control$cov.init=='covOGK')
        Sigma <- covOGK(t(rr), sigmamu = scaleTau2, ...)$wcov
      ## else if (control$cov.init=='QC')
      ##   Sigma <- quadrant.covar(t(rr), ...)$covariance
      else if (control$cov.init=='2SGS') {
        rrNA <- univariate.filter(t(rr), ...)
        Sigma <- GSE(rrNA)@S
      }
      ## } else
      ##   Sigma <- composite.S(t(rr), ...)$S      
    }
    vv <- matrix(0, nrow=p*(p+1)/2, ncol=R)
    lt <- lower.tri(diag(p), diag=TRUE)
    for (r in 1:R) {
      vv[,r] <- as.vector(V[,,r][lt])
    }
    v0 <- as.vector(diag(p)[lt])
    ss <- Sigma[lt]
    if (control$gamma.univ) {
      gamma <- rep(0, R)
      for (r in 1:R)
        gamma[r] <- lm.fit(x=vv[,r,drop=FALSE], y=ss)$coefficients
      if (is.null(eta0))
        eta0 <- lm.fit(x=v0, y=ss)$coefficients
    } else {
      eta <- lm.fit(x=cbind(v0, vv), y=ss)$coefficients
      repeta <- eta < control$epsilon & c(0,control$lower) >=0
      eta[repeta] <- pmax(rep(control$epsilon, R+1), c(0,control$lower))[repeta]
      gamma <- eta[-1]/eta[1]
      if (is.null(eta0))
        eta0 <- eta[1]
    }
    repgamma <- gamma < control$epsilon & control$lower >=0
    gamma[repgamma] <- pmax(rep(control$epsilon, R), control$lower)[repgamma]
  } else {
    if (is.null(eta0)) 
      stop("When you provide 'gamma' you must provide 'eta0' too")
  }
  if (is.null(Sigma))
    Sigma <- eta0*Vprod(V=V, gamma=gamma) ##  V0 is added in Vprod
  
  if (is.null(scales)) {
    RR <- rssr(resid=rr, Sigma=Sigma)
    if (control$psi=="rocke") {
      warning("scales evaluated using optimal rho function with default tuning values")
      
      psi <- "optimal"
      tuning.chi <- 1
      bb <- 0.5
    } else {
      psi <- control$psi
      tuning.chi <- control$tuning.chi
      bb <- control$bb
    }
    scales <- doSsteppw(RR=RR, scale=rep(10, p*(p-1)/2), bb=bb, cc=tuning.chi, psi=psi, tol=control$rel.tol.scale/1000, verbose=(control$trace.lev>2))
  }
  if (is.null(scale)) {  
    Sigmastar <- Sigma/det(Sigma)^(1/nrow(Sigma))
    Sigmastarinv <- solve(Sigmastar)
    RR <- rep(0, n)
    for (i in 1:n)
      RR[i] <- drop(rr[,i]%*%Sigmastarinv%*%rr[,i])
    if (control$psi!="rocke")
      scale <- doSstep(m=RR, scale=10, bb=control$bb, cc=control$tuning.chi, psi=control$psi, tol=control$rel.tol.scale/1000, verbose=(control$trace.lev>2))
    else
      scale <- doSsteprocke(m=RR, scale=10, bb=control$bb, p=p, arp=control$arp.chi, tol=control$rel.tol.scale/1000, verbose=(control$trace.lev>2))
  }

  if (control$trace.lev>1) {
    cat('Initial values \n')
    cat('beta: ', beta, '\n')
    cat('gamma: ', gamma, '\n')
    cat('eta0: ', eta0, '\n\n')
  }
##END# Initial values

  ans <- list(beta=beta, gamma=gamma, eta0=eta0, scales=scales, scale=scale)
  return(ans)
}
