###
### greedyGLM.R
###

# Greedy algorithm to find approximate posterior mode
# Searches for model M that maximizes log-likelihood(theta_M) + log-prior(theta_M|M) + log-prior(M),
# where theta_M: MLE of theta under model M.
# The maximization starts at the null model and sequentially adds/drops any covariate increasing the target function.
# The algorigthm stops when no covariates are added/dropped after a full pass, i.e. a local mode has been found
#
# Note: the algorithm only considers models with at most length(y) variables.
#
# Input:
# - y, x, xadj: response, covariates and adjustment covariates
# - family: family argument passed on to glm
# - priorCoef: prior on regression coefficients
# - priorDelta: prior on model space
# - maxit: maximum number of iterations
# Output: logical vector with length ncol(x) indicating the variables included in the model

greedyGLM <- function(y, x, xadj, family, priorCoef, priorDelta, maxit=100) {
  pluginJoint <- function(sel) {
    p <- sum(sel)
    if (p>0 & p<=length(y)) {
      glm1 <- glm(y ~ x[,sel,drop=FALSE] + xadj -1, family=family)
      ans <- -.5*glm1$deviance + cfprior(matrix(coef(glm1)[1:p],nrow=1)) + modelprior(sel)
    } else if (p==0) {
      glm1 <- glm(y ~ xadj -1, family=family)
      ans <- -.5*glm1$deviance + modelprior(sel)
    } else { ans <- -Inf }
    return(ans)
  }
  #Set prior on model space
  if (priorDelta@priorDistr=='uniform') {
    modelprior <- unifPrior
  } else if (priorDelta@priorDistr=='binomial') {
    if (all(c('alpha.p','beta.p') %in% names(priorDelta@priorPars))) {
      modelprior= function(sel) bbPrior(sel,alpha=priorDelta@priorPars['alpha.p'],beta=priorDelta@priorPars['beta.p'],logscale=TRUE)
    } else modelprior= function(sel) binomPrior(sel, prob=priorDelta@priorPars['p'], logscale=TRUE)
  } else {
    stop('Prior on model space prDelta not recognized')
  }
  #Set coefficients prior
  if (priorCoef@priorDistr=='pMOM') {
    if (all(c('a.tau','b.tau') %in% names(priorCoef@priorPars))) {
      cfprior <- function(th) dmom(th,a.tau=priorCoef@priorPars['a.tau'],b.tau=priorCoef@priorPars['b.tau'],r=priorCoef@priorPars['r'],logscale=TRUE)
    } else {
      cfprior <- function(th) dmom(th,tau=priorCoef@priorPars['tau'],r=priorCoef@priorPars['r'],logscale=TRUE)
    }
  } else if (priorCoef@priorDistr=='piMOM') {
    cfprior <- function(th) dimom(th,tau=priorCoef@priorPars['tau'],logscale=TRUE)
  } else if (priorCoef@priorDistr=='peMOM') {
    if (all(c('a.tau','b.tau') %in% names(priorCoef@priorPars))) {
      cfprior <- function(th) demom(th,a.tau=priorCoef@priorPars['a.tau'],b.tau=priorCoef@priorPars['b.tau'],logscale=TRUE)
    } else {
      cfprior <- function(th) demom(th,tau=priorCoef@priorPars['tau'],logscale=TRUE)
    }
  } else {
    stop('prior on coefficients not recognized')
  }
  #Greedy iterations
  sel <- rep(FALSE,ncol(x))
  mcur <- pluginJoint(sel)
  nchanges <- 1; niter <- 1
  while (nchanges>0 & niter<maxit) {
    nchanges <- 0; niter <- niter+1
    for (i in 1:ncol(x)) {
      selnew <- sel; selnew[i] <- !selnew[i]
      mnew <- pluginJoint(selnew)
      if (mnew>mcur) { sel[i] <- selnew[i]; mcur <- mnew; nchanges <- nchanges+1 }
    }
  }
  return(sel)
}

