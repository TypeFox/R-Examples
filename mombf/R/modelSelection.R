###
### modelSelection.R
###

### Methods for msfit objects

setMethod("show", signature(object='msfit'), function(object) {
  cat('msfit object with',object$p,'variables and',object$family,'residual distribution\n')
  ifelse(any(object$postMode!=0), paste('  Posterior mode: covariate',which(object$postMode==1)), '  Posterior mode: null model')
  cat("Use postProb() to get posterior model probabilities\n")
  cat("Elements $margpp, $postMode, $postSample and $coef contain further information (see help('msfit') and help('modelSelection') for details)\n")
}
)

setMethod("postProb", signature(object='msfit'), function(object, nmax, method='norm') {
  if (method=='norm') {
    modelpp <- unique(data.frame(object$postSample==1, logpp=object$postProb))
    modelpp <- data.frame(modelid= apply(modelpp[,1:(ncol(modelpp)-1)], 1, function(z) paste(which(z),collapse=',')), logpp=modelpp$logpp)
    modelpp$logpp <- modelpp$logpp - modelpp$logpp[1]
    modelpp$pp <- exp(modelpp$logpp)/sum(exp(modelpp$logpp))
  } else if (method=='exact') {
    modelpp <- apply(object$postSample==1, 1, function(z) paste(which(z),collapse=','))
    modelpp <- table(modelpp)/length(modelpp)
    modelpp <- data.frame(modelid=names(modelpp), pp=as.numeric(modelpp))
  } else {
    stop("Argument 'method' not recognized")
  }
  modelpp <- modelpp[order(modelpp$pp,decreasing=TRUE),]
  if (!missing(nmax)) modelpp <- modelpp[1:nmax,]
  if (object$family=='auto') {
    modelid <- as.character(modelpp[,'modelid'])
    twopiece <- laplace <- logical(nrow(modelpp))
    twopiece[grep(as.character(object$p+1),modelid)] <- TRUE
    laplace[grep(as.character(object$p+2),modelid)] <- TRUE
    family <- character(nrow(modelpp))
    family[(!twopiece) & (!laplace)] <- 'normal'
    family[twopiece & (!laplace)] <- 'twopiecenormal'
    family[(!twopiece) & laplace] <- 'laplace'
    family[twopiece & laplace] <- 'twopiecelaplace'
    modelid <- sub(paste(',',object$p+1,sep=''),'',modelid)
    modelid <- sub(as.character(object$p+1),'',modelid)  #for null model
    modelid <- sub(paste(',',object$p+2,sep=''),'',modelid)
    modelid <- sub(as.character(object$p+2),'',modelid)  #for null model
    modelpp <- data.frame(modelid=modelid,family=family,pp=modelpp[,'pp'])
  } else {
    modelpp <- data.frame(modelid=modelpp[,'modelid'],family=object$family,pp=modelpp[,'pp'])
  }
  modelpp[,c('modelid','family','pp')]
}
)



#### General model selection routines
modelSelection <- function(y, x, center=TRUE, scale=TRUE, niter=10^4, thinning=1, burnin=round(niter/10), family='normal', priorCoef=momprior(tau=0.348), priorDelta=modelbbprior(alpha.p=1,beta.p=1), priorVar=igprior(alpha=.01,lambda=.01), priorSkew=momprior(tau=0.348), phi, deltaini=rep(FALSE,ncol(x)), initSearch='greedy', method='auto', B=10^5, verbose=TRUE) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
# - center: if center==TRUE y and x are centered to have zero mean, therefore eliminating the need to include an intercept term in x.
# - scale: if scale==TRUE y and columns in x are scaled to have standard deviation 1
# - niter: number of Gibbs sampling iterations
# - thinning: MCMC thinning factor, i.e. only one out of each thinning iterations are reported. Defaults to thinning=1, i.e. no thinning
# - burnin: number of burn-in MCMC iterations. Defaults to 10% of niter. Set to 0 for no burn-in.
# - family: assumed residual distribution ('normal','twopiecenormal','laplace','twopiecelaplace')
# - priorCoef: prior distribution for the coefficients. Must be object of class 'msPriorSpec' with slot priorType set to 'coefficients'. Possible values for slot priorDistr are 'pMOM', 'piMOM' and 'peMOM'.
# - priorDelta: prior on model indicator space. Must be object of class 'msPriorSpec' with slot priorType set to 'modelIndicator'. Possible values for slot priorDistr are 'uniform' and 'binomial'
# - priorVar: prior on residual variance. Must be object of class 'msPriorSpec' with slot priorType set to 'nuisancePars'. Slot priorDistr must be equal to 'invgamma'.
# - priorSkew: prior on residual skewness parameter. Ignored unless family=='twopiecenormal' or 'twopiecelaplace'
# - phi: residual variance. Typically this is unknown and therefore left missing. If specified argument priorVar is ignored.
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# - initSearch: algorithm to refine deltaini. initSearch=='greedy' uses a greedy Gibbs sampling search. initSearch=='SCAD' sets deltaini to the non-zero elements in a SCAD fit with cross-validated regularization parameter. initSearch=='none' leaves deltaini unmodified.
# - method: method to compute marginal densities. method=='Laplace' for Laplace approx, method=='MC' for Importance Sampling, method=='Hybrid' for Hybrid Laplace-IS (the latter method is only used for piMOM prior with unknown residual variance phi), method='plugin'
# - B: number of samples to use in Importance Sampling scheme. Ignored if method=='Laplace'.
# - verbose: set verbose==TRUE to print iteration progress
# Output: list
# - postSample: posterior samples
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)

  #Check input
  if (!is.vector(y)) { y <- as.double(as.vector(y)) } else { y <- as.double(y) }
  if (!is.matrix(x)) x <- as.matrix(x)
  ct <- (colMeans(x^2)-colMeans(x)^2)==0
  y <- scale(y,center=center,scale=scale); x[,!ct] <- scale(x[,!ct],center=center,scale=scale)
  if (missing(phi)) { knownphi <- as.integer(0); phi <- double(0) } else { knownphi <- as.integer(1); phi <- as.double(phi) }
  p <- ncol(x); n <- length(y)
  if (missing(deltaini)) {
    deltaini <- integer(0); ndeltaini= as.integer(0)
  } else {
    if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
    if (!is.logical(deltaini)) { stop('deltaini must be of type logical') } else { ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(which(deltaini)-1) }
  }
  if (nrow(x)!=length(y)) stop('nrow(x) must be equal to length(y)')
  if (method=='Laplace') {
    method <- as.integer(0)
  } else if (method=='MC') {
    method <- as.integer(1)
  } else if (method=='Hybrid') {
    if ((priorCoef@priorDistr!='piMOM') | (knownphi==1)) {
      warning("method=='Hybrid' is only available for 'piMOM' priors with unknown phi. Using method=='Laplace' instead")
      method <- as.integer(0)
    } else {
      method <- as.integer(2)
    }
  } else if (method=='auto') {
    if (priorCoef@priorDistr!='pMOM') {
      method <- as.integer(0)
    } else {
      method <- as.integer(2)
    }
  } else if (method=='plugin') {
    method <- as.integer(2)
  } else {
    stop("Invalid 'method'")
  }

  #Format arguments for .Call
  niter <- as.integer(niter); burnin <- as.integer(burnin); thinning <- as.integer(thinning); B <- as.integer(B)
  sumy2 <- as.double(sum(y^2)); XtX <- t(x) %*% x; ytX <- as.vector(matrix(y,nrow=1) %*% x)

  if (priorCoef@priorDistr=='pMOM') {
    r <- as.integer(priorCoef@priorPars['r']); prior <- as.integer(0)
  } else if (priorCoef@priorDistr=='piMOM') {
    r <- as.integer(1); prior <- as.integer(1)
  } else if (priorCoef@priorDistr=='peMOM') {
    r <- as.integer(1); prior <- as.integer(2)
    stop('eMOM prior not currently implemented. Try function emomLM instead')
  } else if (priorCoef@priorDistr=='zellner') {
    r <- as.integer(1); prior <- as.integer(3)
  } else {
    stop('Prior specified in priorDistr not recognized')
  }
  tau <- as.double(priorCoef@priorPars['tau'])
  alpha <- as.double(priorVar@priorPars['alpha']); lambda <- as.double(priorVar@priorPars['lambda'])

  taualpha <- as.double(priorSkew@priorPars['tau'])
  if (family=='auto') { family <- 0 } else if (family=='normal') { family <- 1 } else if (family=='twopiecenormal') { family <- 2 } else if (family=='laplace') { family <- 3 } else if (family=='twopiecelaplace') { family <- 4 } else stop("family not available")
  family <- as.integer(family)

  if (priorDelta@priorDistr=='uniform') {
    prDelta <- as.integer(0)
    prDeltap <- as.double(0)
    parprDeltap <- double(2)
  } else if (priorDelta@priorDistr=='binomial') {
    if ('p' %in% priorDelta@priorPars) {
      prDelta <- as.integer(1)
      prDeltap <- as.double(priorDelta@priorPars['p'])
      if ((prDeltap<=0) | (prDeltap>=1)) stop("p must be between 0 and 1 for priorDelta@priorDistr=='binomial'")
      parprDeltap <- double(2)
    } else {
      prDelta <- as.integer(2)
      prDeltap <- as.double(.5)
      parprDeltap <- as.double(priorDelta@priorPars[c('alpha.p','beta.p')])
    }
  } else {
    stop('Prior specified in priorDelta not recognized')
  }


  #Initialize
  if (family==0) { postMode <- rep(as.integer(0),p+1) } else { postMode <- rep(as.integer(0),p) }
  postModeProb <- double(1)
  if (initSearch=='greedy') {
    niterGreed <- as.integer(100)
    if (family==0) { famgreedy <- as.integer(1) } else { famgreedy <- family }
    ans <- .Call("greedyVarSelCI", postMode,postModeProb,knownphi,famgreedy,prior,niterGreed,ndeltaini,deltaini,n,p,y,sumy2,x,XtX,ytX,method,B,alpha,lambda,phi,tau,taualpha,r,prDelta,prDeltap,parprDeltap,as.integer(verbose))
    ndeltaini <- as.integer(sum(postMode)); deltaini <- as.integer(which(as.logical(postMode))-1)
  } else if (initSearch=='SCAD') {
    #require(ncvreg)
    if (verbose) cat("Initializing via SCAD cross-validation...")
    deltaini <- rep(TRUE,ncol(x))
    cvscad <- cv.ncvreg(X=x[,!ct],y=y-mean(y),family="gaussian",penalty="SCAD",nfolds=10,dfmax=1000,max.iter=10^4)
    deltaini[!ct] <- ncvreg(X=x[,!ct],y=y-mean(y),penalty='SCAD',dfmax=1000,lambda=rep(cvscad$lambda[cvscad$cv],2))$beta[-1,1]!=0
    ndeltaini <- as.integer(sum(deltaini)); deltaini <- as.integer(which(deltaini)-1)
    if (verbose) cat(" Done\n")
  }

  #Run MCMC
  mcmc2save <- floor((niter-burnin)/thinning)
  if (family != 0) {
    postSample <- rep(as.integer(0),p*mcmc2save)
    margpp <- double(p)
  } else {
    postSample <- rep(as.integer(0),(p+1)*mcmc2save)
    margpp <- double(p+1)
  }
  if (prDelta==2) postOther <- double(mcmc2save) else postOther <- double(0)
  postProb <- double(mcmc2save)
  ans <- .Call("modelSelectionCI", postSample,postOther,margpp,postMode,postModeProb,postProb,knownphi,family,prior,niter,thinning,burnin,ndeltaini,deltaini,n,p,y,sumy2,as.double(x),XtX,ytX,method,B,alpha,lambda,phi,tau,taualpha,r,prDelta,prDeltap,parprDeltap,as.integer(verbose))
  postSample <- matrix(postSample,ncol=ifelse(family!=0,p,p+1))
  if (family==0) { family <- 'auto' } else if (family==1) { family <- 'normal' } else if (family==2) { family <- 'twopiecenormal' } else if (family==3) { family <- 'laplace' } else if (family==4) { family <- 'twopiecelaplace' }
  if (family=='normal') {
    coef <- rep(0,ncol(x))
    if (any(postMode)) {
      pm <- postMode(y=y,x=x[,postMode[1:ncol(x)]==1,drop=FALSE],priorCoef=priorCoef)
      coef[postMode==1] <- pm$coef
    }
  } else {
    coef <- rep(NA,ncol(x))
  }

  ans <- list(postSample=postSample,postOther=postOther,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb,coef=coef,family=family,p=ncol(x))
  new("msfit",ans)
}

greedymodelSelectionR <- function(y, x, niter=100, marginalFunction, priorFunction, betaBinPrior, deltaini=rep(FALSE,ncol(x)), verbose=TRUE, ...) {
  #Greedy version of modelSelectionR where variables with prob>0.5 at current iteration are included deterministically (prob<.5 excluded)
  p <- ncol(x)
  if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
  if (!missing(betaBinPrior)) {
    #Initialize probBin
    if ((betaBinPrior['alpha.p']>1) && (betaBinPrior['beta.p']>1)) {
      probBin <- (betaBinPrior['alpha.p']-1)/(betaBinPrior['alpha.p']+betaBinPrior['beta.p']-2)
    } else {
      probBin <- (betaBinPrior['alpha.p'])/(betaBinPrior['alpha.p']+betaBinPrior['beta.p'])
    }
    postOther <- matrix(NA,nrow=niter,ncol=1); colnames(postOther) <- c('probBin')
    priorFunction <- function(sel, logscale=TRUE) dbinom(x=sum(sel),size=length(sel),prob=probBin,log=logscale)
  } else {
    postOther <- matrix(NA,nrow=niter,ncol=0)
  }
  #Greedy iterations
  sel <- deltaini
  mcur <- marginalFunction(y=y,x=x[,sel,drop=FALSE],logscale=TRUE,...) + priorFunction(sel,logscale=TRUE)
  nchanges <- 1; itcur <- 1
  nn <- names(x)
  #browser()
  while (nchanges>0 & itcur<niter) {
    nchanges <- 0; itcur <- itcur+1
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      mnew <- marginalFunction(y=y,x=x[,selnew,drop=FALSE],logscale=TRUE,...) + priorFunction(selnew,logscale=TRUE)
      if (mnew>mcur) { sel[i]=selnew[i]; mcur=mnew; nchanges=nchanges+1; if (verbose) cat(paste(ifelse(sel[i],"Added","Dropped"),nn[i],"\n",collapse=" ")) }
    }
  }
  return(sel)
}

modelSelectionR <- function(y, x, niter=10^4, marginalFunction, priorFunction, betaBinPrior, deltaini=rep(FALSE,ncol(x)), verbose=TRUE, ...) {
# Input
# - y: vector with response variable
# - x: design matrix with all potential predictors
# - niter: number of Gibbs sampling iterations
# - marginalFunction: function to compute the marginal density of the data under each model
# - priorFunction: function to compute the model prior probability
# - betaBinPrior: if specified, priorFunction argument is ignored and set to a binomial prior with Beta-hyperprior for the success prob. betaBinPrior should be a vector with named elements 'alpha.p' and 'beta.p', e.g. betaBinPrior= c(alpha.p=1,beta.p=1)
# - deltaini: logical vector of length ncol(x) indicating which coefficients should be initialized to be non-zero. Defaults to all variables being excluded from the model
# ...: other arguments to be passed on to marginalFunction
# Output: list
# - postSample: posterior samples for model indicator
# - postOther: posterior samples for other parameters (probBin: success probability for Binomial prior on number of coefficients in the model)
# - margpp: marginal posterior probability for inclusion of each covariate (approx by averaging marginal post prob for inclusion in each Gibbs iteration. This approx is more accurate than simply taking colMeans(postSample))
# - postMode: model with highest posterior probability amongst all those visited
# - postModeProb: unnormalized posterior prob of posterior mode (log scale)
# - postProb: unnormalized posterior prob of each visited model (log scale)
if (class(y)=='Surv') {
  if ((length(y)/2) != nrow(x)) stop("Dimensions of y and x do not match")
} else {
  if (length(y) != nrow(x)) stop("Dimensions of y and x do not match")
}
if (any(is.na(x))) stop("x cannot have missing values")
p <- ncol(x)
if (length(deltaini)!=p) stop('deltaini must be of length ncol(x)')
if (!missing(betaBinPrior)) {
  #Initialize probBin
  if ((betaBinPrior['alpha.p']>1) && (betaBinPrior['beta.p']>1)) {
    probBin <- (betaBinPrior['alpha.p']-1)/(betaBinPrior['alpha.p']+betaBinPrior['beta.p']-2)
  } else {
    probBin <- (betaBinPrior['alpha.p'])/(betaBinPrior['alpha.p']+betaBinPrior['beta.p'])
  }
  postOther <- matrix(NA,nrow=niter,ncol=1); colnames(postOther) <- c('probBin')
  priorFunction <- function(sel, logscale=TRUE) dbinom(x=sum(sel),size=length(sel),prob=probBin,log=logscale)
} else {
  postOther <- matrix(NA,nrow=niter,ncol=0)
}
if (ncol(x)>12) {
  sel <- postMode <- deltaini
  currentJ <- postModeProb <- marginalFunction(y=y,x=x[,sel,drop=FALSE],logscale=TRUE,...) + priorFunction(sel,logscale=TRUE)
  postSample <- matrix(NA,nrow=niter,ncol=p)
  margpp <- double(p)
  postProb <- double(niter)
  k <- 1; postProb[k] <- postModeProb
  names(postProb)[k] <- paste("V",which(sel),collapse=',',sep='')
  niter10 <- ceiling(niter/10)
  for (i in 1:niter) {
    for (j in 1:p) {
      selnew <- sel; selnew[j] <- !sel[j]
      namenew <- paste(which(selnew),collapse=',')
      newJ <- postProb[namenew]
      if (is.na(newJ)) {
        newJ <- marginalFunction(y=y,x=x[,selnew,drop=FALSE],logscale=TRUE,...) + priorFunction(selnew,logscale=TRUE)
        k <- k+1; postProb[k] <- newJ
        names(postProb)[k] <- paste("V",which(selnew),collapse=',',sep='')
      }
      if (newJ>postModeProb) {
        postModeProb <- newJ
        postMode <- selnew
      }
      pp <- 1/(1+exp(-currentJ+newJ))
      if (sel[j]) {  #if variable in the model
        sel[j] <- runif(1)<pp
        margpp[j] <- margpp[j]+pp
      } else {       #if variable not in the model
        sel[j] <- runif(1)>pp
        margpp[j] <- margpp[j]+1-pp
      }
      if (sel[j]==selnew[j]) {  #if value was updated, update marginal and prior densities
        currentJ <- newJ
      }
    }
    if (!missing(betaBinPrior)) {
      probBin <- rbeta(1,betaBinPrior['alpha.p']+sum(sel), betaBinPrior['beta.p']+sum(!sel))
      postOther[i,'probBin'] <- probBin
    }
    postSample[i,] <- sel
    postProb[i] <- currentJ
    if (verbose && ((i%%niter10)==0)) cat('.')
  }
  margpp <- margpp/niter
  if (verbose) cat('\n')
  #Format postProb
  modelid <- sapply(apply(postSample,1,which),function(z) paste("V",z,collapse=',',sep=''))
  postProb <- postProb[modelid]
} else {
  if (verbose) cat(paste("Computing posterior probabilities for all",2^ncol(x),"models..."))
  models <- expand.grid(lapply(1:ncol(x),function(z) c(FALSE,TRUE)))
  postProb <- apply(models,1, function(z) marginalFunction(y=y,x=x[,z,drop=FALSE],logscale=TRUE,...) + priorFunction(z,logscale=TRUE))
  if (verbose) cat('\n')
  postMode <- models[which.max(postProb),]
  postModeProb <- max(postProb)
  pp <- postProb - postModeProb; pp <- exp(pp)/sum(exp(pp))
  sampledmodels <- rep(1:nrow(models), rmultinom(1,size=niter,prob=pp)[,1])
  postSample <- models[sampledmodels,]
  postProb <- postProb[sampledmodels]
  margpp <- as.vector(t(models) %*% matrix(pp,ncol=1))
}
ans <- list(postSample=postSample,postOther=postOther,margpp=margpp,postMode=postMode,postModeProb=postModeProb,postProb=postProb)
ans <- new("msfit",ans)
return(ans)
}


#Gibbs model selection using BIC to approximate marginal likelihood
# - y, x, xadj: response, covariates under selection and adjustment covariates
# - family: glm family, passed on to glm
# - niter: number of Gibbs iteration
# - burnin: burn-in iterations
# - modelPrior: function evaluating model log-prior probability. Takes a logical vector as input
# Returns:
# - postModel, postCoef1, postCoef2, margpp (analogous to pmomPM. postCoef1 & postCoef2 are MLEs under each visited model)
modelselBIC <- function(y, x, xadj, family, niter=1000, burnin= round(.1*niter), modelPrior, verbose=TRUE) {
  pluginJoint <- function(sel) {
    p <- sum(sel)
    ans <- vector("list",2); names(ans) <- c('marginal','coef')
    if (p>0 & p<=length(y)) {
      glm1 <- glm(y ~ x[,sel,drop=FALSE] + xadj -1, family=family)
      ans$marginal <- -.5*glm1$deviance - .5*log(length(y))*(glm1$df.null-glm1$df.residual) + modelPrior(sel)
      ans$coef1 <- coef(glm1)[1:p]; ans$coef2 <- coef(glm1)[-1:-p]
    } else if (p==0) {
      glm1 <- glm(y ~ xadj -1, family=family)
      ans$marginal <- -.5*glm1$deviance - .5*log(length(y))*(glm1$df.null-glm1$df.residual) + modelPrior(sel)
      ans$coef1 <- double(0); ans$coef2 <- coef(glm1)
    } else { ans$marginal <- -Inf; ans$coef1 <- ans$coef2 <- 0}
    return(ans)
  }
  #Greedy iterations
  if (verbose) cat("Initializing...")
  sel <- rep(FALSE,ncol(x))
  mcur <- pluginJoint(sel)$marginal
  nchanges <- 1; it <- 1
  while (nchanges>0 & it<100) {
    nchanges <- 0; it <- it+1
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      mnew <- pluginJoint(selnew)$marginal
      if (mnew>mcur) { sel[i] <- selnew[i]; mcur <- mnew; nchanges <- nchanges+1 }
    }
  }
  if (verbose) { cat(" Done\nGibbs sampling") }
  #Gibbs iterations
  niter10 <- ceiling(niter/10)
  postModel <- matrix(NA,nrow=niter,ncol=ncol(x))
  postCoef1 <- matrix(0,nrow=niter,ncol=ncol(x))
  postCoef2 <- matrix(0,nrow=niter,ncol=ncol(xadj))
  curmod <- pluginJoint(sel)
  for (j in 1:niter) {
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      newmod <- pluginJoint(selnew)
      mnew <- newmod$marginal
      if (runif(1) < 1/(1+exp(mcur-mnew))) { sel[i] <- selnew[i]; mcur <- mnew; curmod <- newmod }
    }
    postModel[j,] <- sel
    postCoef1[j,sel] <- curmod$coef1; postCoef2[j,] <- curmod$coef2
    if (verbose & ((j%%niter10)==0)) cat(".")
  }
  if (verbose) cat("Done\n")
  #Return output
  if (burnin>0) { postModel <- postModel[-1:-burnin,,drop=FALSE]; postCoef1 <- postCoef1[-1:-burnin,,drop=FALSE]; postCoef2 <- postCoef2[-1:-burnin,,drop=FALSE] }
  ans <- list(postModel=postModel, postCoef1=postCoef1, postCoef2=postCoef2, margpp=colMeans(postModel))
}


## Common prior distributions on model space

binomPrior <- function(sel, prob=.5, logscale=TRUE) {  dbinom(x=sum(sel),size=length(sel),prob=prob,log=logscale) }
unifPrior <- function(sel, logscale=TRUE) { ifelse(logscale,-length(sel)*log(2),2^(-length(sel)))  }
bbPrior <- function(sel, alpha=1, beta=1, logscale=TRUE) {
  ans <- lbeta(sum(sel) + alpha, sum(!sel) + beta) - lbeta(alpha,beta)
  ifelse(logscale,ans,exp(ans))
}


bbPriorTrunc <- function (sel, logscale=TRUE, maxvars=10) {
  #Same as bbPrior with prob=0 when variables > maxvars
  if (sum(sel)<=maxvars) { ans <- bbPrior(sel, logscale=logscale) } else { ans <- ifelse(logscale, -Inf, 0) }
  return(ans)
}

unifPriorTrunc <- function (sel, logscale=TRUE, maxvars=10) {
  #Same as unifPrior with prob=0 when variables > maxvars
  if (sum(sel)<=maxvars) { ans <- unifPrior(sel, logscale=logscale) } else { ans <- ifelse(logscale, -Inf, 0) }
  return(ans)
}
