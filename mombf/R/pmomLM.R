###
### pmomLM.R
###

pmomLM <- function(y, x, xadj, center=FALSE, scale=FALSE, niter=10^4, thinning=1, burnin=round(niter/10), priorCoef, priorDelta, priorVar, initSearch='greedy', verbose=TRUE) {
  #Check input
  if (!is.vector(y)) { y <- as.double(as.vector(y)) } else { y <- as.double(y) }
  if (!is.matrix(x)) x <- as.matrix(x)
  y <- scale(y,center=center,scale=scale)
  ct <- (colMeans(x^2)-colMeans(x)^2)==0; x[,!ct] <- scale(x[,!ct],center=center,scale=scale)
  if (missing(priorCoef)) { priorCoef <- new("msPriorSpec",priorType='coefficients',priorDistr='pMOM',priorPars=c(a.tau=1,b.tau=.135,r=1)) }
  if (missing(priorDelta)) { priorDelta <- new("msPriorSpec",priorType='modelIndicator',priorDistr='uniform',priorPars=double(0)) }
  if (missing(priorVar)) { priorVar <- new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=.01,lambda=.01)) }
  p1 <- ncol(x); n <- length(y)
  if (nrow(x)!=length(y)) stop('nrow(x) must be equal to length(y)')
  if (!missing(xadj)) {
    if (!is.matrix(xadj)) xadj <- as.matrix(xadj)
    p2 <- ncol(xadj)
    ctadj <- (colMeans(xadj^2)-colMeans(xadj)^2)==0; xadj[,!ctadj] <- scale(xadj[,!ctadj],center=center,scale=scale)
    if (nrow(xadj)!=length(y)) stop('nrow(xadj) must be equal to length(y)')
  } else {
    p2 <- as.integer(0); xadj <- double(1)
  }
  
  #Format arguments for .Call
  niter <- as.integer(niter); burnin <- as.integer(burnin); thinning <- as.integer(thinning)
  isbinary <- as.integer(0); ybinary <- integer(0)
  sumy2 <- as.double(sum(y^2)); XtX <- t(x) %*% x; ytX <- as.vector(matrix(y,nrow=1) %*% x); colsumx1sq <- as.double(colSums(x^2))
  if (priorCoef@priorDistr=='pMOM') {
    r <- as.integer(priorCoef@priorPars['r']); prCoef <- as.integer(1)
    if (all(c('a.tau','b.tau') %in% names(priorCoef@priorPars))) {
      atau1 <- as.double(priorCoef@priorPars['a.tau']); btau1 <- as.double(priorCoef@priorPars['b.tau']); tau1 <- as.double(.2); priorTau1 <- as.integer(1)
    } else {
      atau1 <- btau1 <- as.double(0); tau1 <- as.double(priorCoef@priorPars['tau']); priorTau1 <- as.integer(0)
    }
    if (is.na(priorCoef@priorPars['tau.adj'])) tau2 <- as.double(10^6) else tau2 <- priorCoef@priorPars['tau.adj']
  } else {
    stop("When calling pmomLM priorCoef@priorDistr must be equal to 'pMOM'")
  }

  alpha <- as.double(priorVar@priorPars['alpha']); lambda <- as.double(priorVar@priorPars['lambda']) 
  if (priorDelta@priorDistr=='uniform') {
    priorModel <- as.integer(0)
    prModelpar <- as.double(0)
  } else if (priorDelta@priorDistr=='binomial') {
    if ('p' %in% priorDelta@priorPars) {
      priorModel <- as.integer(1)
      prModelpar <- as.double(priorDelta@priorPars['p'])
      if ((prModelpar<=0) | (prModelpar>=1)) stop("p must be between 0 and 1 for priorDelta@priorDistr=='binomial'")
    } else {
      priorModel <- as.integer(2)
      prModelpar <- as.double(priorDelta@priorPars[c('alpha.p','beta.p')])
    }
  } else {
    stop('Prior specified in priorDelta not recognized')
  }

  if (p2>0) {
    S2 <- t(xadj) %*% xadj + diag(1/tau2,nrow=p2); cholS2 <- chol(S2, pivot = TRUE); cholS2 <- cholS2[,order(attr(cholS2, "pivot"))]
    S2inv <- solve(S2); cholS2inv <- chol(S2inv, pivot = TRUE); cholS2inv <- cholS2inv[,order(attr(cholS2inv, "pivot"))]
  } else {
    S2 <- cholS2 <- S2inv <- cholS2inv <- double(1)
  }
  
  #Initialize
  if (initSearch=='greedy') {
    niterGreed <- as.integer(100)
    msfit <- modelSelection(y=y,x=x,center=center,scale=scale,niter=1,priorCoef=priorCoef,priorDelta=priorDelta,priorVar=priorVar,initSearch="greedy",method='Laplace',verbose=FALSE) 
    ndeltaini <- as.integer(sum(msfit$postMode)); deltaini <- as.integer(msfit$postMode)
  } else if (initSearch=='SCAD') {
    #require(ncvreg)
    if (verbose) cat("Initializing via SCAD cross-validation...")
    deltaini <- rep(TRUE,ncol(x))
    cvscad <- cv.ncvreg(X=x[,!ct],y=y-mean(y),family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
    deltaini[!ct] <- ncvreg(X=x[,!ct],y=y-mean(y),penalty='SCAD',dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1,1]!=0
    ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(deltaini)
    if (verbose) cat(" Done\n")
  } else if (initSearch=='none') {
    ndeltaini <- as.integer(0); deltaini <- as.integer(rep(0,ncol(x)))
  }
  if (((ndeltaini+p2)<n) & (p2>0)) {
    iniCoef1 <- rep(0,p1)
    if (ndeltaini>0) {
      lmini <- lm(y ~ -1 + x[,deltaini==1] + xadj)
      iniCoef1[deltaini==1] <- coef(lmini)[1:ndeltaini];
      iniCoef2 <- coef(lmini)[-1:-ndeltaini]  
    } else { lmini <- lm(y ~ -1+xadj); iniCoef2 <- coef(lmini) }
    iniCoef1 <- as.double(iniCoef1)
    iniPhi <- as.double(summary(lmini)$sigma^2)
  } else {
    lmini <- lm(y ~ -1+xadj); iniCoef2 <- coef(lmini) 
    iniCoef1 <- as.double(rep(0,p1))
    iniPhi <- as.double(summary(lmini)$sigma^2)
  }
  if (priorTau1!=0) iniOthers <- atau1/btau1 else iniOthers <- tau1

  #Run MCMC
  mcmc2save <- floor((niter-burnin)/thinning)
  postModel <- integer(p1*mcmc2save); margpp <- double(p1); postCoef1 <- double(p1*mcmc2save); postCoef2 <- double(p2*mcmc2save); postPhi <- double(mcmc2save)
  if (priorTau1!=0) postOther <- double(mcmc2save) else postOther <- double(1)
  x <- as.double(x); xadj <- as.double(xadj)

  ans <- .Call("pmomLM_I", postModel,margpp,postCoef1,postCoef2,postPhi,postOther,niter,thinning,burnin,ndeltaini,deltaini,iniCoef1,iniCoef2,iniPhi,iniOthers,as.integer(verbose),n,p1,p2,isbinary,ybinary,y,sumy2,x,xadj,XtX,ytX,cholS2,S2inv,cholS2inv,colsumx1sq,alpha,lambda,prCoef,r,tau1,tau2,priorTau1,atau1,btau1,priorModel,prModelpar)
  postModel <- matrix(postModel,ncol=p1); postCoef1 <- matrix(postCoef1,ncol=p1)
  if (p2>0) postCoef2 <- matrix(postCoef2,ncol=p2) else postCoef2 <- NA
  if (priorTau1==0) postOther <- NA
  return(list(postModel=postModel,postCoef1=postCoef1,postCoef2=postCoef2,postPhi=postPhi,postOther=postOther,margpp=margpp))
}




pmomLMR <- function(y, x, xadj, phi, r=1, tau, tau.adj=10^6, alpha.phi=.01, lambda.phi=.01, a.tau=1, b.tau=.135, niter=10^3, modelPrior=bbPrior, initSearch='SCAD', verbose=TRUE) {
#Fit linear model with pmom prior on regression coefficients
# Input
# - y: vector with response variable
# - x: design matrix with covariates to be selected
# - xadj: design matrix with ajustment covariates for which no selection process is to be performed (i.e. always included in the model). xadj should include a column of 1's to account for the intercept term. By default xadj is set to matrix(1,ncol=1,nrow=length(y))
# - phi: residual variance. If unknown leave phi missing, a prior phi ~ Inverse Gamma (alpha.phi/2,lambda.phi/2) is used.
# - r: power parameter. pMOM prior for non-zero coefficients is proportional to theta^(2*r) N(theta;0,tau*phi)
# - tau: prior dispersion for pmom prior on the coefficients associated to x
# - tau.adj: prior dispersion for multivariate normal prior on the coefficients associated to xadj
# - alpha.phi, lambda.phi: prior on phi is IG(alpha.phi/2,lambda.phi/2)
# - a.tau, b.tau: if tau unspecified, prior on tau is IG(a.tau/2,b.tau/2). Defaults to values giving 5% prob to interval (-.2,.2)
# - niter: number of Gibbs sampling iterations
# - modelPrior: function to compute the model log-prior probability
# Output: list with 2 elements
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample)).
#require(mvtnorm)
if (missing(phi)) { unknownPhi <- TRUE } else { unknownPhi <- FALSE }
if (missing(tau)) { unknownTau <- TRUE } else { unknownTau <- FALSE }
if (is.vector(y)) y <- matrix(y,ncol=1)
if (missing(xadj)) xadj <- matrix(1,nrow=nrow(y),ncol=1)
#Pre-compute useful quantities
n <- nrow(y); p1 <- ncol(x); p2 <- ncol(xadj)
XtX <- t(x) %*% x
S2 <- t(xadj) %*% xadj + diag(1/tau.adj,nrow=p2)
S2inv <- solve(S2)
cholS2inv <- chol(S2inv, pivot = TRUE)
cholS2inv <- cholS2inv[,order(attr(cholS2inv, "pivot"))]
#Initialize
postDelta <- postTheta1 <- matrix(NA,nrow=niter,ncol=p1)
postTheta2 <- matrix(NA,nrow=niter,ncol=p2)
postPhi <- postTau <- double(niter)
if (initSearch=='none') {
  sel <- rep(FALSE,p1)
  postDelta[1,] <- sel
  postTheta1[1,] <- rep(0,p1)
} else if (initSearch=='SCAD') {
  #require(ncvreg)
  cvscad <- cv.ncvreg(X=x,y=y,family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
  postTheta1[1,] <- ncvreg(X=x,y=y,penalty="SCAD",dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1, 1]
  postDelta[1,] <- postTheta1[1,]!=0
}
ndeltaini <- sum(postDelta[1,])
if (((ndeltaini+p2)<n) & (p2>0)) {
  if (ndeltaini>0) {
    lmini <- lm(y ~ -1 + x[,postDelta[1,]] + xadj) 
    postTheta1[1,postDelta[1,]] <- coef(lmini)[1:ndeltaini]
    postTheta2[1,] <- coef(lmini)[-1:-ndeltaini]
  } else { lmini <- lm(y ~ -1 + xadj); postTheta1[1,] <- rep(0,p1); postTheta2[1,] <- coef(lmini) }
} else {
  if (ndeltaini>0) {
    lmini <- lm(y ~ -1 + x[,postDelta[1,]])
    postTheta1[1,postDelta[1,]] <- coef(lmini)
    postTheta2[1,] <- double(p2)
  } else { lmini <- lm(y ~ -1 + xadj); postTheta1[1,] <- rep(0,p1); postTheta2[1,] <- coef(lmini) }
}
linpred1 <- x %*% t(postTheta1[1,,drop=FALSE])
linpred2 <- xadj %*% t(postTheta2[1,,drop=FALSE])
e <- y-linpred1-linpred2
postPhi[1] <- ifelse(unknownPhi, summary(lmini)$sigma^2, phi)
postTau[1] <- ifelse(unknownTau, .2, tau)
#Iterate
for (i in 2:niter) {
  #Sample delta1, theta1
  curDelta <- postDelta[i-1,]; curTheta1 <- postTheta1[i-1,]
  for (j in 1:p1) {
    ej <- e+curTheta1[j]*x[,j]
    newval <- MHTheta1pmom(ej,j=j,delta=curDelta,theta1=curTheta1,phi=postPhi[i-1],r=r,tau=postTau[i-1],xj=x[,j],padj=p2,modelPrior=modelPrior)
    curDelta[j] <- newval$delta; curTheta1[j] <- newval$theta1
    if (newval$accept) e <- ej - curTheta1[j]*x[,j]   #Update residuals
  }
  postDelta[i,] <- curDelta; postTheta1[i,] <- curTheta1
  #Sample theta2
  e <- e+linpred2
  postTheta2[i,] <- simTheta2(e=e, xadj=xadj, S2inv=S2inv, cholS2inv=cholS2inv, phi=postPhi[i-1])
  linpred2 <- xadj %*% t(postTheta2[i,,drop=FALSE])
  e <- e - linpred2
  #Sample phi
  postPhi[i] <- ifelse(unknownPhi, simPhipmom(alpha.phi=alpha.phi,lambda.phi=lambda.phi,n=n,r=r,delta=postDelta[i,],p2=p2,theta1=postTheta1[i,],theta2=postTheta2[i,],tau=postTau[i-1],tau.adj=tau.adj,ssr=sum(e^2)), phi)
  #Sample tau
  postTau[i] <- ifelse(unknownTau, simTaupmom(a.tau=a.tau,b.tau=b.tau,r=r,delta=postDelta[i,],theta1=postTheta1[i,],phi=postPhi[i]), tau)
  if (verbose & ((i %% (niter/10))==0)) cat('.')
}
if (verbose) cat('\n')
ans <- cbind(postDelta,postTheta1,postTheta2,postPhi,postTau)
colnames(ans) <- c(paste('delta',1:ncol(postDelta),sep=''),paste('theta',1:ncol(postTheta1),sep=''),paste('thetaAdj',1:ncol(postTheta2),sep=''),'phi','tau')
return(ans)
}


## Routines for the R implementation of the scheme

proposalpmom <- function(m,S,phi,r,tau,e,xj,m1,nu) {
#Approximate univariate pmom posterior with a 2 component T mixture with nu df
# Posterior: N(e; xj*theta; phi*I) * pmom(theta; phi, r, tau) / m1
# Mixture: w1 * T_nu(theta;mu1,sigma21) + (1-w1) * T_nu(theta;mu2,sigma22)  (sigma21, sigma22 denote variances)
# - m,S: posterior parameters
# - phi: residual variance
# - r: pmom prior power is 2*r
# - tau: prior dispersion parameter
# - e: response variable
# - xj: predictor
# - m1: normalization constant
# - nu: desired degrees of freedom
# Output: named vector with parameters of approximating mixture
  mu <- .5*(m + c(-1,1)*sqrt(m^2+8*r*phi/S)) #Posterior mode
  fmode <- exp(sum(dnorm(e, mu[2]*xj, sd=sqrt(phi), log=TRUE)) + dmom(mu[2],tau=tau,phi=phi,r=r,logscale=TRUE) - m1) #Value at mode
  sigma2 <- 1/diag(fppmomNeg(mu,m=m,S=S,phi=phi,tau=tau,r=r)) #Proposal variances
  ct2 <- exp(lgamma(.5*nu+.5) - .5*log(nu) - lgamma(.5*nu) - .5*log(pi*sigma2[2]))
  w1 <- max(0,(fmode - ct2)/(dnorm(mu[2],mu[1],sd=sqrt(sigma2[1])) - ct2)) #Weight
  ans <- c(mu,sigma2,w1)
  names(ans) <- c('mu1','mu2','sigma21','sigma22','w1')
  return(ans)
}

MHTheta1pmom <- function(e,j,delta,theta1,phi,r,tau,xj,padj,modelPrior) {
  #MH step to simulate (delta[j], theta1[j]) from its posterior given the data, delta[-j], theta1[-j], theta2 and phi parameters
  # Input
  # - e: partial residuals, i.e. y - predicted y given all covariates except covariate j
  # - j: index of the element in delta and theta1 to update
  # - delta: current value for delta
  # - theta1: current value for theta1
  # - phi: current value for phi
  # - r: pmom prior power is 2*r
  # - tau: current value for tau
  # - xj: vector containing the column of the design matrix associated with theta1[j]
  # - padj: number of adjustment covariates. Models with total number of vars >= never accepted
  # - modelPrior: function to compute the model log-prior probability
  # Ouput: list with the following elements
  # - delta: new value for delta[j] (can be the same as input value if proposal not accepted)
  # - theta1: new value for theta1[j]
  # - accept: logical variable indicated whether proposed new value has been accepted or not
  #Propose
  pcur <- sum(delta)
  m1 <- pmomMargKuniv(y=e, x=xj, phi=phi, tau=tau, logscale=TRUE)
  logbf <- sum(dnorm(e,0,sd=sqrt(phi),log=TRUE)) - m1
  delta0 <- delta1 <- delta; delta0[j] <- FALSE; delta1[j] <- TRUE
  logpratio01 <- modelPrior(delta0) - modelPrior(delta1)
  if ((delta[j]==0) & ((pcur+padj) >= length(e))) {
      p <- 0
  } else {
      p <- 1/(1 + exp(logbf+logpratio01))
  }
  deltaProp <- rbinom(n=1,size=1,prob=p)
  nu <- floor(sqrt(ifelse(is.matrix(e),nrow(e),length(e))))
  #Acceptance prob
  if ((!delta[j]) & (deltaProp==0)) {
    thetaProp <- 0
    lambda <- 1
  } else {
    S <- sum(xj^2) + 1/tau; m <- sum(xj*e)/S
    propPars <- proposalpmom(m=m,S=S,phi=phi,r=r,tau=tau,e=e,xj=xj,m1=m1,nu=nu)
    thetaProp <- rTmix2comp(pars=propPars, df=nu)
    if (delta[j] & (deltaProp==1)) {
      lhood <- sum(dnorm(e,thetaProp*xj,sd=sqrt(phi),log=TRUE)) - sum(dnorm(e,theta1[j]*xj,sd=sqrt(phi),log=TRUE))
      lprior <- dmom(thetaProp,tau=tau,phi=phi,r=r,logscale=TRUE) - dmom(theta1[j],tau=tau,phi=phi,r=r,logscale=TRUE)
      lprop <- dTmix2comp(theta1[j],pars=propPars,df=nu,logscale=TRUE) - dTmix2comp(thetaProp,pars=propPars,df=nu,logscale=TRUE)
      lambda <- exp(lhood + lprior + lprop)
    } else if ((!delta[j]) & (deltaProp==1)) {
      num <- sum(dnorm(e,thetaProp*xj,sd=sqrt(phi),log=TRUE)) + dmom(thetaProp,tau=tau,phi=phi,r=r,logscale=TRUE)
      den <- dTmix2comp(thetaProp,pars=propPars,df=nu,logscale=TRUE) + m1
      lambda <- exp(num-den)
    } else if ((delta[j]) & (deltaProp==0)) {
      thetaProp <- 0
      num <- dTmix2comp(theta1[j],pars=propPars,df=nu,logscale=TRUE) + m1
      den <- sum(dnorm(e,theta1[j]*xj,sd=sqrt(phi),log=TRUE)) + dmom(theta1[j],tau=tau,phi=phi,r=r,logscale=TRUE)
      lambda <- exp(num-den)
    }
  }
  if (runif(1)<lambda) {
    ans <- list(delta=(deltaProp==1), theta1=thetaProp, accept=TRUE)
  } else {
    ans <- list(delta=delta[j], theta1=theta1[j], accept=FALSE)
  }
  return(ans)
}

simPhipmom <- function(alpha.phi, lambda.phi, n, r, delta, p2, theta1, theta2, tau, tau.adj, ssr) {
 #Draw from posterior of the variance given all other parameters under a pmom prior
  a <- alpha.phi + n + (2*r+1)*sum(delta) + p2
  b <- lambda.phi + sum(theta1^2)/tau + sum(theta2^2)/tau.adj + ssr
  1/rgamma(1,a/2,b/2)
}

simTaupmom <- function(a.tau, b.tau, r, delta, theta1, phi) {
  #Draw from posterior of tau given all other parameters under a pmom prior
  a <- a.tau + (2*r+1)*sum(delta)
  b <- b.tau + sum(theta1^2)/phi
  1/rgamma(1,a/2,b/2)
}

simTheta2 <- function(e, xadj, S2inv, cholS2inv, phi) {
  #Simulate Theta2 ~ N(S2inv %*% t(xadj) %*% e, phi*S2inv)
  m2 <- S2inv %*% t(xadj) %*% e
  sweep(matrix(rnorm(ncol(xadj)),nrow=1) %*% cholS2inv * sqrt(phi), 2, m2, "+")
}

dTmix2comp <- function(th, pars, df, logscale=TRUE) {
  #Evaluate density of 2 component T with df degrees of freedom mixture
  ans <- pars['w1']*dmvt(th,pars['mu1'],matrix(pars['sigma21'],nrow=1),df=df,log=FALSE) + (1-pars['w1'])*dmvt(th,pars['mu2'],matrix(pars['sigma22'],nrow=1),df=df,log=FALSE)
  if (logscale) ans <- log(ans)
  return(ans)
}

rTmix2comp <- function(pars, df) {
  #Generate single draw from T mixture with 2 components and df degrees of freedom
  ifelse(runif(1)<pars['w1'], pars['mu1']+rmvt(n=1,sigma=matrix(pars['sigma21'],nrow=1),df=df), pars['mu2']+rmvt(n=1,sigma=matrix(pars['sigma22'],nrow=1),df=df))
}

pmomMargKuniv <- function(y,x,phi,tau=1,r=1,logscale=TRUE) {
#Univariate marginal density under a product MOM prior (known variance case)
# integral N(y; x*theta, phi*I) * (theta^2/(tau*phi))^r * N(theta; 0; tau*phi) / (2r-1)!! d theta
# - y: response variable (must be a vector)
# - x: design matrix (must be a vector)
# - phi: residual variance
# - tau: prior variance parameter (defaults to length(y))
# - logscale: if set to TRUE the log of the integral is returned
  #require(actuar)
  n <- length(y)
  if (n != length(y)) stop("Dimensions of x and y don't match")
  if (missing(tau)) tau <- n
  s <- sum(x^2) + 1/tau
  m <- sum(x*y)/s
  I <- log(actuar::mnorm(2*r, mean=m, sd=sqrt(phi/s)))
  ans <- I -.5*(sum(y^2) - s*m^2)/phi - .5*n*log(2*pi*phi) - .5*(log(s)+log(tau)) - sum(log(seq(1,2*r-1,by=2))) - r*log(tau*phi)
  if (!logscale) ans <- exp(ans)
  return(ans)
}

