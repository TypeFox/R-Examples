# Generate component failure times given parameter value and censoring status
# censoring: vector with 0=uncensored, -1=left censored, 1=right censored
rexpCompGivenParm <- function(pars, failTime, censoring, ...) {
  rate <- pars$rate
  y <- vector("numeric", length(censoring))
  for(cens in 1:length(censoring)) {
    if(censoring[cens]==0) {
      y[cens] = failTime
      next
    } else if(censoring[cens]==-1) {
      # solve y=(1 - Exp[-l x])/(1 - Exp[-l T]) for x
      y[cens] = log(1/(runif(1)*(exp(-rate*failTime)-1)+1))/rate
    } else {
      # by Exponential memoryless
      y[cens] = failTime + rexp(1, rate)
    }
  }
  y
}


##### IID #####
# Extra ... parameters
#   priorShape & priorScale

# Random draw from posterior of Exponential given data and Gamma prior
rexpParmGivenDataIID <- function(data, priorShape, priorScale) {
  posteriorShape <- priorShape + length(data)
  posteriorScale <- 1/(1/priorScale + sum(data))
  
  rgamma(1, shape=posteriorShape, scale=posteriorScale)
}

# t         = vector of system/network failure times
# signature = signature of the system/network for which inference is performed
# priorShape = shape parameter of Gamma prior
# priorScale = scale parameter of Gamma prior
maskedInferenceIIDExponential <- function(t, signature, iter, priorShape, priorScale) {
  maskedInferenceIIDCustom(t, signature, function(q, pars, ...){pexp(q, pars$rate)}, function(x, pars, ...){dexp(x, pars$rate)}, rexpParmGivenDataIID, rexpCompGivenParm, c(rate=rgamma(1, shape=priorShape, scale=priorScale)), iter, priorShape=priorShape, priorScale=priorScale)
}


##### EXCHANGEABLE #####
# Extra ... parameters to Log-Normal priors on shape and scale
#   priorMuShape, priorSigmaShape, priorMuScale, priorSigmaScale

# Supporting fns
#library(mcmc, lib.loc="~/R/library")
#library(numDeriv, lib.loc="~/R/library")
alphaBetaLPost <- function(musig, mm, nn, priorMu_Mu, priorSigma_Mu, priorMu_Sigma, priorSigma_Sigma, data) {
  mu <- musig[1]
  sig <- musig[2]
  if(mu <= 0 || sig <= 0) return(-Inf)
  res <- (-nn * (mu^2/sig) * log((sig/mu)) +
    nn * lgamma((mu^2/sig) + mm) -
    nn * lgamma((mu^2/sig)) -
    (priorMu_Mu - log(mu))^2 / (2 * priorSigma_Mu^2) -
    log(mu) -
    (priorMu_Sigma - log(sig))^2 / (2 * priorSigma_Sigma^2) -
    log(sig) -
    ((mu^2/sig) + mm) * sum(log(1/(sig/mu) + rowSums(data))) -
    log(sig))
  #print(c(mu,sig,res))
  res
}
alphaBetaLPost2 <- function(mu, sig, mm, nn, priorMu_Mu, priorSigma_Mu, priorMu_Sigma, priorSigma_Sigma, data) {
  alphaBetaLPost(c(mu, sig), mm, nn, priorMu_Mu, priorSigma_Mu, priorMu_Sigma, priorSigma_Sigma, data)
}

# Random draw from posterior of Exponential given data and Gamma prior
rexpParmGivenDataEXCH <- function(data, priorMu_Mu, priorSigma_Mu, priorMu_Sigma, priorSigma_Sigma) {
  alphaBetaLPost3 <- Vectorize(alphaBetaLPost2, c("mu", "sig"))
  mm <- dim(data)[2]
  nn <- dim(data)[1]
  
  #cat("Data matrix is",nn,"x",mm)
  
  # p(a,b|y)
  #axisM <- seq(1,11,length.out=100)
  #axisS <- seq(0.1,10,length.out=100)
  #z <- outer(axisM,axisS,FUN="alphaBetaLPost3", mm=mm, nn=nn, mu_alpha=priorMuShape, sigma_alpha=priorSigmaShape, mu_beta=priorMuScale, sigma_beta=priorSigmaScale, data=data)
  #contour(axisM,axisS,exp(z-max(z)), nlevels=50)
  
  ###opt <- optim(c(last[nn+1]*last[nn+2],last[nn+1]*last[nn+2]^2), alphaBetaLPost, m=m, nn=nn, mu_alpha=priorMuShape, sigma_alpha=priorSigmaShape, mu_beta=priorMuScale, sigma_beta=priorSigmaScale, x=data, method="L-BFGS-B", lower=c(0.001,0.001), upper=c(10,10), control=list(fnscale=-1), hessian=TRUE)
  ##opt <- optim(c(last[nn+1]*last[nn+2],last[nn+1]*last[nn+2]^2), alphaBetaLPost, mm=mm, nn=nn, mu_alpha=priorMuShape, sigma_alpha=priorSigmaShape, mu_beta=priorMuScale, sigma_beta=priorSigmaScale, data=data, control=list(fnscale=-1), hessian=TRUE)
  #opt <- optim(c(priorMuShape*priorMuScale,priorMuShape*priorMuScale^2), alphaBetaLPost, mm=mm, nn=nn, mu_alpha=priorMuShape, sigma_alpha=priorSigmaShape, mu_beta=priorMuScale, sigma_beta=priorSigmaScale, data=data, control=list(fnscale=-1), hessian=TRUE)
  opt <- optim(c(priorMu_Mu, priorMu_Sigma), alphaBetaLPost, mm=mm, nn=nn, priorMu_Mu=priorMu_Mu, priorSigma_Mu=priorSigma_Mu, priorMu_Sigma=priorMu_Sigma, priorSigma_Sigma=priorSigma_Sigma, data=data, control=list(fnscale=-1), hessian=TRUE)
  #print(opt)
  muhat <- opt$par
  #muhat <- c(last[nn+1]*last[nn+2],last[nn+1]*last[nn+2]^2)
  H <- -opt$hessian
  #H <- -hessian(alphaBetaLPost, muhat, method="Richardson", mm=mm, nn=nn, mu_alpha=priorMuShape, sigma_alpha=priorSigmaShape, mu_beta=priorMuScale, sigma_beta=priorSigmaScale, data=data)
  #print(H)
  res <- metrop(alphaBetaLPost, muhat, 100, scale=chol(solve(H)), mm=mm, nn=nn, priorMu_Mu=priorMu_Mu, priorSigma_Mu=priorSigma_Mu, priorMu_Sigma=priorMu_Sigma, priorSigma_Sigma=priorSigma_Sigma, data=data)
  #points(res$batch[100,1],res$batch[100,2], col="red", pch=20)
  parm <- c(res$batch[100,1]^2/res$batch[100,2], res$batch[100,2]/res$batch[100,1])
  #points(parm[1]*parm[2],parm[1]*parm[2]*parm[2], col="green", pch=20)
  
  # p(psi|a,b,y)
  #print(parm[1]+mm)
  #print(parm[2]+rowSums(data))
  #c(rgamma(nn, shape=parm[1]+mm, scale=1/(1/parm[2]+rowSums(data))), parm)
  list(parm, as.list(rgamma(nn, shape=parm[1]+mm, scale=1/(1/parm[2]+rowSums(data)))))
}

# t         = vector of system/network failure times
# signature = signature of the system/network for which inference is performed
maskedInferenceEXCHExponential <- function(t, signature, iter, priorMu_Mu, priorSigma_Mu, priorMu_Sigma, priorSigma_Sigma) {
  print({priorMu <- exp(priorMu_Mu + priorSigma_Mu*rnorm(1))})
  print({priorSigma <- exp(priorMu_Sigma + priorSigma_Sigma*rnorm(1))})
  print({priorShape <- priorMu^2/priorSigma})
  print({priorScale <- priorSigma/priorMu})
  print({rateInit <- rgamma(1, shape=priorShape, scale=priorScale)})
  maskedInferenceEXCHCustom(t, signature, function(q, pars, ...){pexp(q, pars$rate)}, function(x, pars, ...){dexp(x, pars$rate)}, rexpParmGivenDataEXCH, rexpCompGivenParm, lapply(1:length(t), function(i) { c(rate=rateInit) }), c(priorShape=priorShape, priorScale=priorScale), iter, priorMu_Mu=priorMu_Mu, priorSigma_Mu=priorSigma_Mu, priorMu_Sigma=priorMu_Sigma, priorSigma_Sigma=priorSigma_Sigma)
}
