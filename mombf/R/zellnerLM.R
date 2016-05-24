zellnerPM <- function(y, x, xadj, tau, tau.adj=10^6, a.tau=1, b.tau=.135, niter=10^3, burnin=round(niter/10), modelPrior=bbPrior, initSearch='SCAD', verbose=TRUE) {
#Fit probit model with zellner prior on regression coefficients
# Input
# - y: vector with response variable (must be a factor with 2 levels or a character which can be converted to factor with 2 levels)
# - x: design matrix with covariates to be selected
# - xadj: design matrix with ajustment covariates for which no selection process is to be performed (i.e. always included in the model). xadj should include a column of 1's to account for the intercept term. By default xadj is set to matrix(1,ncol=1,nrow=length(y))
# - tau: prior dispersion for zellner prior on the coefficients associated to x
# - tau.adj: prior dispersion for multivariate normal prior on the coefficients associated to xadj
# - a.tau, b.tau: if tau unspecified, prior on tau is IG(a.tau/2,b.tau/2). Defaults to values giving 5% prob to interval (-.2,.2)
# - niter: number of Gibbs sampling iterations
# - modelPrior: function to compute the model log-prior probability
# Output: list with 2 elements
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample)).
#require(mvtnorm)
if (missing(tau)) { unknownTau <- TRUE; stop("Case with unknown tau is currently not implemented. Please specify tau") } else { unknownTau <- FALSE }
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
    newval <- MHTheta1zellner(ej,j=j,delta=curDelta,theta1=curTheta1,phi=1,tau=postTau[i-1],xj=x[,j],padj=p2,modelPrior=modelPrior)
    curDelta[j] <- newval$delta; curTheta1[j] <- newval$theta1
    if (newval$accept) e <- ej - curTheta1[j]*x[,j]   #Update residuals
  }
  postDelta[i,] <- curDelta; postTheta1[i,] <- curTheta1
  linpred1 <- x %*% matrix(curTheta1,ncol=1)
  #Sample theta2
  e <- e+linpred2
  postTheta2[i,] <- simTheta2(e=e, xadj=xadj, S2inv=S2inv, cholS2inv=cholS2inv, phi=1)
  linpred2 <- xadj %*% t(postTheta2[i,,drop=FALSE])
  #Sample tau
  postTau[i] <- tau
  #postTau[i] <- ifelse(unknownTau, simTauzellner(a.tau=a.tau,b.tau=b.tau,r=r,delta=postDelta[i,],theta1=postTheta1[i,],phi=1), tau)
  if (verbose & ((i %% (niter/10))==0)) cat('.')
}
if (verbose) cat('\n')
ans <- cbind(postDelta,postTheta1,postTheta2,postTau)
colnames(ans) <- c(paste('delta',1:ncol(postDelta),sep=''),paste('theta',1:ncol(postTheta1),sep=''),paste('thetaAdj',1:ncol(postTheta2),sep=''),'tau')
return(ans[-1:-burnin,,drop=FALSE])
}




MHTheta1zellner <- function(e,j,delta,theta1,phi,tau,xj,padj,modelPrior) {
  #MH step to simulate (delta[j], theta1[j]) from its posterior given the data, delta[-j], theta1[-j], theta2 and phi parameters
  # Input
  # - e: partial residuals, i.e. y - predicted y given all covariates except covariate j
  # - j: index of the element in delta and theta1 to update
  # - delta: current value for delta
  # - theta1: current value for theta1
  # - phi: current value for phi
  # - tau: current value for tau
  # - xj: vector containing the column of the design matrix associated with theta1[j]
  # - padj: number of adjustment covariates (forced inclusion in the model)
  # - modelPrior: function to compute the model log-prior probability
  # Ouput: list with the following elements
  # - delta: new value for delta[j] (can be the same as input value if proposal not accepted)
  # - theta1: new value for theta1[j]
  # - accept: logical variable indicated whether proposed new value has been accepted or not
  #Propose
  pcur <- sum(delta)
  m1 <- zellnerMargKuniv(y=e, x=xj, phi=phi, tau=tau, logscale=TRUE)
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
    thetaProp <- rnorm(1,m,sd=1/sqrt(S))
    if (delta[j] & (deltaProp==1)) {
      lhood <- sum(dnorm(e,thetaProp*xj,sd=sqrt(phi),log=TRUE)) - sum(dnorm(e,theta1[j]*xj,sd=sqrt(phi),log=TRUE))
      lprior <- dnorm(thetaProp,0,sd=sqrt(tau*phi),log=TRUE) - dnorm(theta1[j],0,sd=sqrt(tau*phi),log=TRUE)
      lprop <- dnorm(theta1[j],m,sd=1/sqrt(S),log=TRUE) - dnorm(thetaProp,m,sd=1/sqrt(S),log=TRUE)
      lambda <- exp(lhood + lprior + lprop)
    } else if ((!delta[j]) & (deltaProp==1)) {
      num <- sum(dnorm(e,thetaProp*xj,sd=sqrt(phi),log=TRUE)) + dnorm(thetaProp,0,sd=sqrt(tau*phi),log=TRUE)
      den <- dnorm(thetaProp,m,sd=1/sqrt(S),log=TRUE) + m1
      lambda <- exp(num-den)
    } else if ((delta[j]) & (deltaProp==0)) {
      thetaProp <- 0
      num <- dnorm(theta1[j],m,sd=1/sqrt(S),log=TRUE) + m1
      den <- sum(dnorm(e,theta1[j]*xj,sd=sqrt(phi),log=TRUE)) + dnorm(theta1[j],0,sd=sqrt(tau*phi),log=TRUE)
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


zellnerMargKuniv <- function(y,x,phi,tau=1,logscale=TRUE) {
#Univariate marginal density under Zellner's prior (known variance case)
# integral N(y; x*theta, phi*I) * (theta^2/(tau*phi))^r * N(theta; 0; tau*phi) / (2r-1)!! d theta
# - y: response variable (must be a vector)
# - x: design matrix (must be a vector)
# - phi: residual variance
# - tau: prior variance parameter (defaults to length(y))
# - logscale: if set to TRUE the log of the integral is returned
  n <- length(y)
  if (n != length(y)) stop("Dimensions of x and y don't match")
  if (missing(tau)) tau <- n
  s <- sum(x^2) * (1+1/tau)
  m <- sum(x*y)/s
  ans <- -.5*(sum(y^2) - s*m^2)/phi - .5*n*log(2*pi*phi) - .5*(log(s)+log(tau)) 
  if (!logscale) ans <- exp(ans)
  return(ans)
}



####### MARGINAL LIKELIHOOD FOR ZELLNER'S PRIOR

zellnerMarginalK <- function(sel, y, x, phi, tau, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under Zellner's prior
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) x <- matrix(x,ncol=1)
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel)); 
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  phi <- as.double(phi); tau <- as.double(tau)
  logscale <- as.integer(logscale)
  ans <- .Call("zellnerMarginalKI", sel, nsel, n, p, y, sumy2, XtX, ytX, phi, tau, logscale)
  return(ans)
}


zellnerMarginalKR <- function(y, x, phi, tau, logscale=TRUE) {
  #Marginal likelihood for Zellner's prior 
  n <- length(y); p <- ncol(x)
  if (p==0) {
    ans <- sum(dnorm(y,0,sd=sqrt(phi),log=TRUE))
  } else {
    S <- (1+1/tau) * t(x) %*% x
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    ans <- -.5*(sum(y^2) - t(m) %*% S %*% m)/phi - .5*n*log(2*pi*phi) - .5*p*log(tau) - log(sqrt(det(S)))
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}


zellnerMarginalU <- function(sel, y, x, alpha=0.001, lambda=0.001, tau=1, logscale=TRUE, XtX, ytX) {
#Marginal density of the data y~N(x*theta,phi*I) under Zellner's prior (unknown variance)
  if (is.matrix(y)) y <- as.vector(y)
  if (is.vector(x)) { x <- matrix(x,ncol=1) } else { x <- as.matrix(x) }
  if (missing(XtX)) { XtX <- t(x) %*% x } else { XtX <- as.matrix(XtX) }
  if (missing(ytX)) { ytX <- as.vector(matrix(y,nrow=1) %*% x) } else { ytX <- as.vector(ytX) }
  if (is.logical(sel)) sel <- which(sel)
  if ((length(sel)>0) && ((min(sel)<1) | max(sel)>ncol(x))) stop('Invalid specification of parameter sel')
  sel <- as.integer(sel-1); nsel <- as.integer(length(sel)); 
  p <- as.integer(ncol(x)); n <- as.integer(nrow(x))
  y <- as.double(y); sumy2 <- sum(y^2)
  tau <- as.double(tau)
  logscale <- as.integer(logscale)
  alpha <- as.double(alpha); lambda <- as.double(lambda)
  ans <- .Call("zellnerMarginalUI",sel,nsel,n,p,y,sumy2,x,XtX,ytX,tau,logscale,alpha,lambda)
  return(ans);
}

zellnerMarginalUR <- function(y, x, alpha=0.001, lambda=0.001, tau, logscale=TRUE) {
  #Marginal likelihood for Zellner's prior
  n <- length(y); p <- ncol(x)
  if (ncol(x)==0) {
    term1 <- .5*(n + alpha)
    num <- .5*alpha*log(lambda) + lgamma(term1)
    den <- .5*n*log(pi) + lgamma(alpha/2)
    ans <- num -den - term1*log(lambda + sum(y^2))
  } else {
    S <- (1+1/tau) * (t(x) %*% x)
    m <- solve(S) %*% t(x) %*% matrix(y,ncol=1)
    nu <- n + alpha
    ss <- as.numeric(lambda + sum(y^2) - t(m) %*% S %*% m)
    #
    num <- lgamma(nu/2) + .5*alpha*log(lambda/2) + .5*nu*log(2) - .5*nu*log(ss)
    den <- .5*n*log(2*pi) + .5*log(det(S)) + (.5*p)*log(tau) + lgamma(alpha/2)
    ans <- num - den
  }
  if (!logscale) ans <- exp(ans)
  return(ans)
}



################################################################################
## Zellner's Bayes factors for linear model objects
################################################################################

###
### zellnerbf.R
###

zellnerbf <- function(lm1, coef, g, theta0, logbf=FALSE) {
  UseMethod("zellnerbf")
}

###
### zellnerbf.lm.R
###

zellnerbf.lm <- function(lm1, coef, g, theta0, logbf=FALSE) {
  if (missing(g)) { stop("'g' must be specified") }
  if (missing(theta0)) {
    theta0 <- rep(0, length(coef))
  } else if (length(theta0) != length(coef)) {
    stop("'theta0' must have the same length as 'coef'")
  }
  thetahat <- lapply(g, function(gg) coef(lm1) * gg/(gg+1))
  V <- lapply(g, function(gg) summary(lm1)$cov.unscaled * gg/(gg+1))
  n <- length(lm1$residuals); p <- length(thetahat[[1]]); p1 <- length(coef)
  if ((min(coef)<1) | (max(coef)>p)) {
    stop("'coef' values must be between 1 and the number of coefficients in 'lm1'")
  }
  ssr <- sum(residuals(lm1)^2); sr <- sqrt(ssr/(n-p))
  bf.zellner <- mapply(function(tt,VV,gg) zbfunknown(tt[coef],VV[coef,coef],n=n,nuisance.theta=p-p1,g=gg,theta0=theta0,ssr=ssr,logbf=logbf), thetahat,V,g)
  bf.zellner
}


###
### zbfunknown.R
###

zbfunknown <- function(theta1hat, V1, n, nuisance.theta, g=1, theta0, ssr, logbf=FALSE) {
  if (missing(theta0)) { theta0 <- rep(0, length(theta1hat)) }
  p1 <- length(theta1hat); p <- p1 + nuisance.theta
  l <- theta1hat-theta0; l <- matrix(l, nrow=1) %*% solve(V1) %*% matrix(l, ncol=1)  
  sigma2hat <- (ssr + l/(1+n*g))/(n-nuisance.theta)
  bf <- (-(n-nuisance.theta)/2)*log(1+n*g*ssr/(ssr+l)) + ((n-p)/2)*log(1+n*g)
  if (!logbf) { bf <- exp(bf) }
  bf
}


###
### zbfknown.R
###

zbfknown <- function(theta1hat, V1, n, g=1, theta0, sigma, logbf=FALSE) {
  if (missing(sigma)) { stop("'sigma' must be specified") }
  if (missing(theta0)) { theta0 <- rep(0, length(theta1hat)) }
  p1 <- length(theta1hat)
  l <- theta1hat-theta0
  l <- matrix(l, nrow=1) %*% solve(V1) %*% matrix(l, ncol=1) * n*g/((1+n*g)*sigma^2) #noncentr param
  muk <- p1+l
  t1 <- matrix(theta1hat-theta0, nrow=1) %*% solve(V1) %*% matrix(theta1hat-theta0, ncol=1) * n*g/((1+n*g)*sigma^2)
  bf <- .5*t1 - .5*p1*log(1+n*g)
  if (!logbf) { bf <- exp(bf) }
  bf
}


###
### zellnerbf.knownsig.R
###

#Bayes factor based on Zellner's g-prior for linear models (known sigma^2 case).
# - theta1hat: vector with estimated value of the coefficients to be tested
# - V1: submatrix of covariance corresponding to elements in theta1hat
# - n: sample size used to fit the model
# - g: prior parameter
# - theta0: hypothesized value for theta1hat (defaults to 0)
# - sigma: residual standard deviation, assumed to be known
# - logbf: if log==TRUE, the natural logarithm of the Bayes factor is returned
zellnerbf.knownsig <- function(theta1hat,
                               V1,
                               n,
                               g=1,
                               theta0,
                               sigma,
                               logbf=FALSE) {
  if (missing(sigma)) {
    stop("'sigma' must be specified")
  }
  if (missing(theta0)) {
    theta0 <- rep(0, length(theta1hat))
  }
  p1 <- length(theta1hat)
  l <- theta1hat-theta0
  l <- matrix(l, nrow=1) %*% solve(V1) %*% matrix(l, ncol=1) * n*g/((1+n*g)*sigma^2) #noncentr param
  muk <- p1+l
  t1 <- matrix(theta1hat-theta0, nrow=1) %*% solve(V1) %*% matrix(theta1hat-theta0, ncol=1) * n*g/((1+n*g)*sigma^2)
  bf <- .5*t1
  if (!logbf) {
    bf <- exp(bf)
  }
  bf
}


###
### zellnerbf.unknownsig.R
###

#Bayes factor based on Zellner's g-prior for linear models (unknown sigma^2 case). 
# - theta1hat: vector with estimated value of the coefficients to be tested
# - V1: submatrix of covariance corresponding to elements in theta1hat
# - n: sample size used to fit the model
# - g: prior parameter
# - theta0: hypothesized value for theta1hat (defaults to 0)
# - ssr: sum of squared residuals
# - logbf: if log==TRUE, the natural logarithm of the Bayes factor is returned
zellnerbf.unknownsig <- function(theta1hat,
                                 V1,
                                 n,
                                 nuisance.theta,
                                 g=1,
                                 theta0,
                                 ssr,
                                 logbf=FALSE) {
  if (missing(theta0)) {
    theta0 <- rep(0, length(theta1hat))
  }
  p1 <- length(theta1hat)
  p <- p1 + nuisance.theta
  l <- theta1hat-theta0
  l <- matrix(l, nrow=1) %*% solve(V1) %*% matrix(l, ncol=1)  
  sigma2hat <- (ssr + l/(1+n*g)) / (n-nuisance.theta)
  muk <- p1+ l* n*g/((1+n*g)*sigma2hat)
  bf <- (-(n-nuisance.theta)/2)*log(1+n*g*ssr/(ssr+l)) + ((n-p)/2)*log(1+n*g)
  if (!logbf) {
    bf <- exp(bf)
  }
  bf
}

