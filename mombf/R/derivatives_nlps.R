###########################################################################################
##
##  GRADIENT AND HESSIAN FOR LOG NON-LOCAL PRIOR (MOM, EMOM AND IMOM PRIOR) DENSITIES. ALSO FOR JOINT (NLP, INV. GAMMA)
##
##  NOTE: THESE ROUTINES ARE ALSO IMPLEMENTED IN C
##
###########################################################################################


###########################################################################################
## DERIVATIVES OF MOM PRIOR
###########################################################################################

#Gradient of log-pMOM(th;0,phi*tau) density wrt th
dmomgrad <- function(th, logphi, tau) { 2/th - th/(exp(logphi)*tau) }

#Hessian of log-pMOM(th;0,phi*tau) density wrt th
dmomhess <- function(th, logphi, tau) { ans <- diag(length(th)); diag(ans) <- -2/th^2 - 1/(exp(logphi)*tau); return(ans) }

#Gradient of log-pMOM(th;0,phi*tau) + log-IG(th;phi,alpha/2,lambda/2) wrt (th, logphi) where logphi=log(phi)
dmomiggrad <- function(th, logphi, tau, alpha, lambda) {
  if (length(th)>0) {
    ans <- c(dmomgrad(th, logphi=logphi, tau=tau), -1.5*length(th) - 0.5*alpha -1 + 0.5*(sum(th)^2/tau + lambda) * exp(-logphi))
  } else {
    ans <- - 0.5*alpha -1 + 0.5*lambda*exp(-logphi)
  }
  return(ans)
}

#Hessian of log-pMOM(0,phi*tau) + log-IG(phi,alpha,lambda) wrt (th, logphi) where logphi=log(phi)
dmomighess <- function(th, logphi, tau, alpha, lambda) {
  p <- length(th)
  ans <- matrix(NA,nrow=p+1,ncol=p+1)
  if (length(th)>0) {
    ans[1:p,1:p] <- dmomhess(th, logphi, tau)
    ans[1:p,p+1] <- ans[p+1,1:p] <- th/(exp(logphi)*tau)
  }
  ans[p+1,p+1] <- -0.5 * exp(-logphi) * (sum(th^2)/tau+lambda)
  return(ans)
}



###########################################################################################
## DERIVATIVES OF EMOM PRIOR
###########################################################################################

#Gradient of log-peMOM(0,phi*tau) density
demomgrad <- function(th, logphi, tau) { 2*tau*exp(logphi)/th^3 - th*exp(-logphi)/tau }

#Hessian of log-peMOM(0,phi*tau) density
demomhess <- function(th, logphi, tau) { ans <- diag(length(th)); diag(ans) <- -6*tau*exp(logphi)/th^4 - exp(-logphi)/tau; return(ans) }

#Gradient of log-peMOM(0,phi*tau) + log-IG(phi,alpha/2,lambda/2) under log(phi) parameterization
demomiggrad <- function(th, logphi, tau, alpha, lambda) {
  if (length(th)>0) {
    ans <- c(demomgrad(th, logphi=logphi, tau=tau), -0.5*length(th) - 0.5*alpha -1 + 0.5*(sum(th)^2/tau + lambda) * exp(-logphi) - exp(logphi)*tau*sum(1/th^2))
  } else {
    ans <- -0.5*alpha -1 + 0.5*lambda * exp(-logphi)
  }
  return(ans)
}

#Hessian of log-peMOM(0,phi*tau) + log-IG(phi,alpha,lambda) under log(phi) parameterization
demomighess <- function(th, logphi, tau, alpha, lambda) {
  p <- length(th)
  ans <- matrix(NA,nrow=p+1,ncol=p+1)
  if (p>0) {
    ans[1:p,1:p] <- demomhess(th, logphi, tau)
    ans[1:p,p+1] <- ans[p+1,1:p] <- th/(exp(logphi)*tau) + 2*tau*exp(logphi)/th^3
  }
  ans[p+1,p+1] <- -0.5 * exp(-logphi) * (sum(th^2)/tau+lambda) - tau * exp(logphi) * sum(1/th^2)
  return(ans)
}


###########################################################################################
## DERIVATIVES OF IMOM PRIOR
###########################################################################################

#Gradient of log-pimom(0,phi*tau) density
dimomgrad <- function(th, logphi, tau) { 2*tau*exp(logphi)/th^3 - 2/th }

#Hessian of log-pimom(0,phi*tau) density
dimomhess <- function(th, logphi, tau) { ans <- diag(length(th)); diag(ans) <- -6*tau*exp(logphi)/th^4 + 2/th^2; return(ans) }

#Gradient of log-pimom(0,phi*tau) + log-IG(phi,alpha/2,lambda/2) under log(phi) parameterization
dimomiggrad <- function(th, logphi, tau, alpha, lambda) {
  if (length(th)>0) {
    ans <- c(dimomgrad(th, logphi=logphi, tau=tau), 0.5*length(th) - 0.5*alpha -1 + 0.5*lambda*exp(-logphi) - exp(logphi)*tau*sum(1/th^2))
  } else {
    ans <- -0.5*alpha -1 + 0.5*lambda*exp(-logphi)
  }
  return(ans)
}

#Hessian of log-pimom(0,phi*tau) + log-IG(phi,alpha,lambda) under log(phi) parameterization
dimomighess <- function(th, logphi, tau, alpha, lambda) {
  p <- length(th)
  ans <- matrix(NA,nrow=p+1,ncol=p+1)
  if (p>0) {
    ans[1:p,1:p] <- dimomhess(th, logphi, tau)
    ans[1:p,p+1] <- ans[p+1,1:p] <- 2*tau*exp(logphi)/th^3
  }
  ans[p+1,p+1] <- -0.5*exp(-logphi)*lambda - tau * exp(logphi) * sum(1/th^2)
  return(ans)
}
