###
### pmomPM.R
###

pmomPM <- function(y, x, xadj, niter=10^4, thinning=1, burnin=round(niter/10), priorCoef, priorDelta, initSearch='greedy', verbose=TRUE) {
  priorVar <- new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=.01,lambda=.01))
  #Check input
  if (is.logical(y) | is.numeric(y)) {
    y <- as.integer(y)
  } else if (is.factor(y)) {
    warning(paste('Converted factor to numeric: 0=',levels(y)[1],' and 1=',levels(y)[2],sep=''))
    y <- as.integer(y)
    y <- as.integer(y - min(y))
  } else {
    stop('y must be either logical, a factor or numeric')
  }
  if (length(unique(y))!=2) stop('y must have two levels')
  if (!is.matrix(x)) x <- as.matrix(x)
  if (missing(priorCoef)) { priorCoef <- new("msPriorSpec",priorType='coefficients',priorDistr='pMOM',priorPars=c(a.tau=1,b.tau=.135,r=1)) }
  if (missing(priorDelta)) { priorDelta <- new("msPriorSpec",priorType='modelIndicator',priorDistr='uniform',priorPars=double(0)) }
  if (missing(priorVar)) { priorVar <- new("msPriorSpec",priorType='nuisancePars',priorDistr='invgamma',priorPars=c(alpha=.01,lambda=.01)) }
  p1 <- ncol(x); n <- length(y)
  if (nrow(x)!=length(y)) stop('nrow(x) must be equal to length(y)')
  if (!missing(xadj)) {
    if (!is.matrix(xadj)) xadj <- as.matrix(xadj)
    p2 <- ncol(xadj)
    if (nrow(xadj)!=length(y)) stop('nrow(xadj) must be equal to length(y)')
  } else {
    p2 <- as.integer(0); xadj <- double(1)
  }
  if (p1>1000) warning('ncol(x)>1000 may require substantial memory and computation time. If you experience problems, specify a smaller number of covariates.')
  
  #Format arguments for .Call
  niter <- as.integer(niter); burnin <- as.integer(burnin); thinning <- as.integer(thinning)
  isbinary <- as.integer(1); ybinary <- as.integer(y); y <- as.double(y)
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
    if (verbose) cat("Initializing via greedy search...")
    msfit <- greedyGLM(y=y,x=x,xadj=xadj,family=binomial(link='probit'),priorCoef=priorCoef,priorDelta=priorDelta,maxit=50)
    ndeltaini <- as.integer(sum(msfit)); deltaini <- as.integer(msfit)
    if (verbose) cat(" Done\n")
  } else if (initSearch=='SCAD') {
    #require(ncvreg)
    if (verbose) cat("Initializing via SCAD cross-validation...")
    deltaini <- rep(TRUE,ncol(x))
    cvscad <- cv.ncvreg(X=x,y=y-mean(y),family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
    deltaini <- ncvreg(X=x,y=y-mean(y),penalty='SCAD',dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1,1]!=0
    ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(deltaini)
    if (verbose) cat(" Done\n")
  } else if (initSearch=='none') {
    ndeltaini <- as.integer(0); deltaini <- as.integer(rep(0,ncol(x)))
  }
  if (((ndeltaini+p2)<n) & (p2>0)) {
    iniCoef1 <- rep(0,p1)
    if (ndeltaini>0) {
      lmini <- glm(y ~ -1 + x[,deltaini==1] + xadj, family=binomial(link='probit'))
      iniCoef1[deltaini==1] <- coef(lmini)[1:ndeltaini];
      iniCoef2 <- coef(lmini)[-1:-ndeltaini]  
    } else { lmini <- glm(y ~ -1+xadj, binomial(link='probit')); iniCoef2 <- coef(lmini) }
    iniCoef1 <- as.double(iniCoef1)
    iniPhi <- as.double(summary(lmini)$sigma^2)
  } else {
    lmini <- glm(y ~ -1+xadj, binomial(link='probit')); iniCoef2 <- coef(lmini) 
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



pmomPMR <- function(y, x, xadj, r=1, tau, tau.adj=10^6, a.tau=1, b.tau=.135, niter=10^3, modelPrior=bbPrior, initSearch='SCAD', verbose=TRUE) {
#Fit probit model with pmom prior on regression coefficients
# Input
# - y: vector with response variable (must be a factor with 2 levels or a character which can be converted to factor with 2 levels)
# - x: design matrix with covariates to be selected
# - xadj: design matrix with ajustment covariates for which no selection process is to be performed (i.e. always included in the model). xadj should include a column of 1's to account for the intercept term. By default xadj is set to matrix(1,ncol=1,nrow=length(y))
# - r: power parameter. pMOM prior for non-zero coefficients is proportional to theta^(2*r) N(theta;0,tau*phi)
# - tau: prior dispersion for pmom prior on the coefficients associated to x
# - tau.adj: prior dispersion for multivariate normal prior on the coefficients associated to xadj
# - a.tau, b.tau: if tau unspecified, prior on tau is IG(a.tau/2,b.tau/2). Defaults to values giving 5% prob to interval (-.2,.2)
# - niter: number of Gibbs sampling iterations
# - modelPrior: function to compute the model log-prior probability
# Output: list with 2 elements
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample)).
#require(mvtnorm)
if (missing(tau)) { unknownTau <- TRUE } else { unknownTau <- FALSE }
if (is.character(y)) { y <- as.numeric(factor(y))-1 } else if (is.factor(y)) { y <- as.numeric(y)-1 }
if (length(unique(y))>2) stop('y has more than 2 levels')
if (missing(xadj)) xadj <- matrix(1,nrow=nrow(y),ncol=1)
#Pre-compute useful quantities
n <- length(y); p1 <- ncol(x); p2 <- ncol(xadj)
XtX <- t(x) %*% x
S2 <- t(xadj) %*% xadj + diag(1/tau.adj,nrow=p2)
S2inv <- solve(S2)
cholS2inv <- chol(S2inv, pivot = TRUE)
cholS2inv <- cholS2inv[,order(attr(cholS2inv, "pivot"))]
#Initialize
postDelta <- postTheta1 <- matrix(NA,nrow=niter,ncol=p1)
postTheta2 <- matrix(NA,nrow=niter,ncol=p2)
if (initSearch=='none') {
  if (verbose) cat("Initializing to null model\n")
  sel <- rep(FALSE,p1)
  postDelta[1,] <- sel
  postTheta1[1,] <- rep(0,p1)
} else if (initSearch=='SCAD') {
  if (verbose) cat("Initializing via SCAD cross-validation")
  #require(ncvreg)
  warn <- options('warn')$warn; options(warn= -1)
  cvscad <- cv.ncvreg(X=x,y=y,family="binomial",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
  postTheta1[1,] <- ncvreg(X=x,y=y,penalty="SCAD",dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1, 1]
  postDelta[1,] <- postTheta1[1,]!=0
  options(warn=warn)
}
postTheta2[1,] <- coef(glm(factor(y) ~ -1 + xadj, family=binomial(link='probit')))
linpred1 <- x %*% t(postTheta1[1,,drop=FALSE])
linpred2 <- xadj %*% t(postTheta2[1,,drop=FALSE])
postTau <- double(niter); postTau[1] <- ifelse(unknownTau, .2, tau)
#Iterate
if (verbose) cat("\nRunning MCMC")
for (i in 2:niter) {
  #Sample latent variables z
  linpred <- linpred1+linpred2; plinpred <- pnorm(-linpred)
  u <- ifelse(y,runif(n,plinpred,1),runif(n,0,plinpred))
  e <- qnorm(u)
  z <- linpred + e
  #Sample delta1, theta1
  curDelta <- postDelta[i-1,]; curTheta1 <- postTheta1[i-1,]
  for (j in 1:p1) {
    ej <- e+curTheta1[j]*x[,j]
    newval <- MHTheta1pmom(ej,j=j,delta=curDelta,theta1=curTheta1,phi=1,r=r,tau=postTau[i-1],xj=x[,j],padj=p2,modelPrior=modelPrior)
    curDelta[j] <- newval$delta; curTheta1[j] <- newval$theta1
    if (newval$accept) e <- ej - curTheta1[j]*x[,j]   #Update residuals
  }
  postDelta[i,] <- curDelta; postTheta1[i,] <- curTheta1
  linpred1 <- x %*% matrix(curTheta1,ncol=1)
  #Sample theta2
  e <- e+linpred2
  postTheta2[i,] <- simTheta2(e=e, xadj=xadj, S2inv=S2inv, cholS2inv=cholS2inv, phi=1)
  linpred2 <- xadj %*% t(postTheta2[i,,drop=FALSE])
  #e <- e - linpred2
  #Sample tau
  postTau[i] <- ifelse(unknownTau, simTaupmom(a.tau=a.tau,b.tau=b.tau,r=r,delta=postDelta[i,],theta1=postTheta1[i,],phi=1), tau)
  if (verbose & ((i %% (niter/10))==0)) cat('.')
}
if (verbose) cat('\n')
ans <- cbind(postDelta,postTheta1,postTheta2,postTau)
colnames(ans) <- c(paste('delta',1:ncol(postDelta),sep=''),paste('theta',1:ncol(postTheta1),sep=''),paste('thetaAdj',1:ncol(postTheta2),sep=''),'tau')
return(ans)
}

