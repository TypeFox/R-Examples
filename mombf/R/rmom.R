##################################################################################
## Routines to simulate from MOM prior and posterior
##################################################################################

setMethod("rnlp", signature(y='ANY',x='matrix',m='missing',V='missing',msfit='msfit'), function(y, x, m, V,msfit, priorCoef, priorVar=igprior(alpha=0.01,lambda=0.01), niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
  if (msfit$family != 'normal') stop("Posterior sampling only implemented for Normal residuals")
  if (!(class(y) %in% c('numeric','Surv'))) stop("y must be of class 'numeric' or 'Surv'")
  #Draw model
  pp <- postProb(msfit,method=pp)
  modelid <- strsplit(as.character(pp$modelid), split=',')
  ndraws <- as.numeric(rmultinom(1, size=niter, prob=pp$pp))
  sel <- ndraws>0; modelid <- modelid[sel]; ndraws <- ndraws[sel]
  #Draw coefficients
  idx <- c(0,cumsum(ndraws))
  if (class(y) == 'numeric') { ##Linear model
    ans <- matrix(0, nrow=niter, ncol=ncol(x)+1)
    for (i in 1:length(modelid)) {
      b <- min(50, ceiling((burnin/niter) * ndraws[i]))
      colsel <- as.numeric(modelid[[i]])
      ans[(idx[i]+1):idx[i+1],c(colsel,ncol(ans))] <- rnlp(y=y, x=x[,colsel,drop=FALSE], priorCoef=priorCoef, priorVar=priorVar, niter=ndraws[i]+b, burnin=b)
    }
    if (is.null(colnames(x))) colnames(ans) <- c(paste('beta',1:ncol(x),sep=''),'phi') else colnames(ans) <- c(colnames(x),'phi')
  } else if (class(y) == 'Surv') {
    ans <- matrix(0, nrow=niter, ncol=ncol(x))
    for (i in 1:length(modelid)) {
      b <- min(50, ceiling((burnin/niter) * ndraws[i]))
      colsel <- as.numeric(modelid[[i]])
      if (length(colsel)>0) {
        ans[(idx[i]+1):idx[i+1],colsel] <- rnlp(y=y, x=x[,colsel,drop=FALSE], priorCoef=priorCoef, priorVar=priorVar, niter=ndraws[i]+b, burnin=b)
      }
    }
    if (is.null(colnames(x))) colnames(ans) <- paste('beta',1:ncol(x),sep='') else colnames(ans) <- colnames(x)
  }
  return(ans)
}
)


setMethod("rnlp", signature(y='ANY',x='matrix',m='missing',V='missing',msfit='missing'), function(y, x, m, V, msfit, priorCoef, priorVar=igprior(alpha=0.01,lambda=0.01), niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
  tau <- as.double(priorCoef@priorPars['tau'])
  if (class(y) == 'numeric') {  ##Linear model
    p <- ncol(x); n <- length(y)
    if (nrow(x) != n) stop('Dimensions of y and x do not match')
    if (priorVar@priorDistr=='invgamma') {
        a_phi <- as.double(priorVar@priorPars['alpha'])
        b_phi <- as.double(priorVar@priorPars['lambda'])
    } else stop("Only invgamma prior for residual variance is currently implemented")
    if (priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) {
      if (priorCoef@priorDistr=='pMOM') {
        prior <- as.integer(0); r <- as.integer(priorCoef@priorPars['r'])
      } else if (priorCoef@priorDistr=='piMOM') {
        prior <- as.integer(1); r <- as.integer(0)
      } else {
        prior <- as.integer(2); r <- as.integer(0)
      }
      if (p==0) {
        ans <- matrix(1/rgamma((niter-burnin)/thinning, .5*(a_phi+n), .5*(b_phi+sum(y^2))), ncol=1)
        colnames(ans) <- 'phi'
      } else {
        ans <- .Call("rnlpPostCI_lm",as.integer(niter),as.integer(burnin),as.integer(thinning),as.double(y),as.double(x),as.integer(p),as.integer(r),tau,a_phi,b_phi,prior)
        ans <- matrix(ans,ncol=p+1)
        if (is.null(colnames(x))) colnames(ans) <- c(paste('beta',1:ncol(x),sep=''),'phi') else colnames(ans) <- c(colnames(x),'phi')
      }
    } else if (priorCoef@priorDistr=='zellner') {
      if (p==0) {
        ans <- matrix(1/rgamma((niter-burnin)/thinning, .5*(a_phi+n), .5*(b_phi+sum(y^2))), ncol=1)
        colnames(ans) <- 'phi'
      } else {
        S <- solve((1+1/tau) * t(x) %*% x)
        m <- as.vector(S %*% t(x) %*% matrix(y,ncol=1))
        ssr <- sum(y * (y - (x %*% m)))
        phi <- 1 / rgamma((niter - burnin)/thinning, 0.5*(n+a_phi), 0.5*(ssr + b_phi))
        beta <- matrix(rnorm(ncol(x)*(niter-burnin)/thinning),ncol=ncol(x)) %*% t(chol(S)) * phi
        beta <- t(t(beta)+m)
        ans <- cbind(beta, phi)
        if (is.null(colnames(x))) colnames(ans) <- c(paste('beta',1:ncol(x),sep=''),'phi') else colnames(ans) <- c(colnames(x),'phi')
      }
    } else stop("This kind of prior is not implemented")
  } else if (class(y)=='Surv') {  #Cox model
    p <- ncol(x)
    if (p==0) {
      ans <- matrix(double(0),nrow=(niter-burnin)/thinning,ncol=0)
    } else {
      fit <- coxph(y ~ ., data=data.frame(x))
      thhat <- matrix(coef(fit),ncol=1); Vinv <- solve(fit$var)
      if (priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) {
        Sinv <- Vinv + diag(p)/tau
      } else if (priorCoef@priorDistr=='zellner') {
        Sinv <- Vinv * (1+1/tau)
      } else stop("This kind of prior is not implemented")
      S <- solve(Sinv)
      m <- S %*% Vinv %*% thhat
      ans <- rnlp(m=as.vector(m), V=S, priorCoef=priorCoef, niter=niter, burnin=burnin, thinning=thinning)
      if (is.null(colnames(x))) colnames(ans) <- paste('beta',1:ncol(x),sep='') else colnames(ans) <- colnames(x)
    }
  } else {
    stop("y must be of class 'numeric' or 'Surv'")
  }
  return(ans)
}
)


setMethod("rnlp", signature(y='missing',x='missing',m='numeric',V='matrix',msfit='missing'), function(y, x, m, V, msfit, priorCoef, priorVar=igprior(alpha=0.01,lambda=0.01), niter=10^3, burnin=round(niter/10), thinning=1, pp='norm') {
  p <- ncol(V)
  tau <- as.double(priorCoef@priorPars['tau'])
  if (priorCoef@priorDistr %in% c('pMOM','peMOM','piMOM')) {
    if (priorCoef@priorDistr=='pMOM') {
      prior <- as.integer(0); r <- as.integer(priorCoef@priorPars['r'])
    } else if (priorCoef@priorDistr=='piMOM') {
      prior <- as.integer(1); r <- as.integer(0)
    } else if (priorCoef@priorDistr=='peMOM') {
      prior <- as.integer(2); r <- as.integer(0)
    } else stop("This kind of prior is not implemented")
    ans <- .Call("rnlpCI",as.integer(niter),as.integer(burnin),as.integer(thinning),as.double(m),as.double(V),as.integer(p),as.integer(r),tau,prior)
    ans <- matrix(ans,ncol=p)
  } else if (priorCoef@priorDistr == 'zellner') {
    ans <- matrix(rnorm(p*(niter-burnin)/thinning),nrow=p) %*% t(chol(V))
  } else stop("This kind of prior is not implemented")
  if (is.null(names(m))) colnames(ans) <- paste('beta',1:length(m),sep='') else colnames(ans) <- names(m)
  return(ans)
}
)



##################################################################################
## Routines to simulate from truncated Normal
##################################################################################

rtnorm <- function(n, m, sd, lower, upper) {
  if (length(lower) != length(upper)) stop('length of lower and upper must match')
  if (length(lower)>0) {
    lower[1] <- max(-1.0e10,lower[1])
    upper[length(upper)] <- min(1.0e10,upper[length(upper)])
    ans <- .Call("rnorm_truncMultCI",as.integer(n),as.double(lower),as.double(upper),as.double(m),as.double(sd))
  } else {
    ans <- rnorm(n, m, sd)
  }
  return(ans)
}

rtmvnorm <- function(n, m, Sigma, SigmaInv, lower, upper, within, method='Gibbs', burnin=round(.1*n)) {
# Multivariate normal samples under rectangular constraint
# Input
# - n: number of draws
# - m: multivariate normal mean
# - Sigma: multivariate normal covariance
# - lower: vector with lower truncation points
# - upper: vector with upper truncation points
# - within: if TRUE, each variable is truncated to be >=lower and <= upper. If FALSE, it's truncated to be <lower or >upper
# - method: set method=='Gibbs' for Gibbs sampling, and method=='MH' for independent proposal MH
# - burnin: number of burn-in iterations
# Output: n draws obtained via Gibbs sampling after orthogonalization
  if (length(lower)==1) lower <- rep(1,length(m))
  if (length(upper)==1) upper <- rep(1,length(m))
  if (length(lower)!=length(m)) stop('Length of lower and m do not match')
  if (length(upper)!=length(m)) stop('Length of upper and m do not match')
  if (nrow(Sigma)!=length(m) | ncol(Sigma)!=length(m)) stop('Dimensions of m and Sigma do no match')
  if (!(method %in% c('Gibbs','MH'))) stop('Method should be Gibbs or MH')
  method <- as.integer(ifelse(method=='Gibbs',1,2))
  ans <- .Call("rtmvnormCI",as.integer(n), as.double(m), as.double(Sigma), as.double(lower), as.double(upper), as.integer(within), method)
  matrix(ans,ncol=length(m))
}

rtmvnormProd <- function(n, m, Sigma, k=1, lower=0, upper=Inf, burnin=round(.1*n)) {
# Multivariate normal samples under product contraint  lower <= prod(x^k) <= upper
# Input
# - m: multivariate normal mean
# - Sigma: multivariate normal covariance
# - k: power of product
# - lower: lower truncation point
# - upper: upper truncation point
# - burnin: number of burn-in iterations
# Output: n draws obtained via Gibbs sampling
  lowtrunc <- ifelse(lower==0,as.integer(0),as.integer(1))
  uptrunc <- ifelse(upper==Inf,as.integer(0),as.integer(1))
  if (upper==Inf) { uptrunc <- as.integer(0); upper <- as.double(0) } else { uptrunc <- as.integer(1); upper <- as.double(upper) }
  ans= .Call("rtmvnormProdCI",as.integer(n), as.double(m), as.double(Sigma), as.integer(k), as.double(lower), upper, lowtrunc, uptrunc, as.integer(burnin));
  matrix(ans,nrow=n)
}

#rmvnormTrunc0 <- function(n, m, S, Sinv, t0) {
#  if (length(m) != nrow(S)) stop('Dimensions of m and S must match')
#  if (length(t0)==1) t0 <- as.double(rep(t0,length(m)))
#  if (length(t0)!= length(m)) stop('Dimensions of m and t0 must match')
#  if (missing(Sinv)) Sinv <- solve(S)
#  ans <- .Call("rmvnorm_trunc0R", as.integer(n), as.double(m), as.double(Sinv), as.double(t0))
#  matrix(ans,nrow=n)
#}
