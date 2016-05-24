# evaluating h,H,S as well as alpha (scale parameter in survreg), lp (linear predictor)
bsgw.eval <- function(coeff, X, Xs, t, ordweib=FALSE, alpha.min=0.1, alpha.max=10.0) {
  K <- ncol(X)
  Ks <- ncol(Xs)
  beta <- coeff[1:K] # scale coefficients, corresponds to ordinary Weibull coefficients
  betas <- coeff[K+1:Ks] # shape coefficients
  lp <- X%*%beta
  if (ordweib) {
    alpha <- exp(betas[1])
  } else {
    alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-Xs%*%betas))
  }
  expterm <- exp(lp)
  H <- t^(alpha) * expterm
  h <- alpha * t^(alpha-1) * expterm
  S <- exp(-H)
  return (list(h=h, H=H, S=S, scale=1/alpha, lp=lp))
}

bsgw.loglike <- function(coeff, X, Xs, t, s, ...) {
  ret <- bsgw.eval(coeff, X, Xs, t, ...)
  return (sum((s*log(ret$h)-ret$H)))
}

bsgw.logpost <- function(coeff, X, Xs, t, s, lambda, lambdas, ...) {
  return (bsgw.loglike(coeff, X, Xs, t, s, ...) - lambda*sum(abs(coeff[1:ncol(X)])) - lambdas*sum(abs(coeff[ncol(X)+1:ncol(Xs)])))
}

bsgw.mcmc <- function(X, Xs, t, s, lambda, lambdas, iter=1000, sd.thresh=1e-06
                      , init=NULL, ordweib=FALSE, alpha.min=0.1, alpha.max=10.0, beta.max=log(20), betas.max=5.0
                      , print.level=2, nskip=10) {
  mbeta <- 100 # TODO: convert from magic number to parameter that can be controlled
  w <- 1.0 # TODO: same as above
  K <- ncol(X)
  Ks <- ncol(Xs)
  beta_smp <- array(NA, dim=c(iter,K))
  betas_smp <- array(NA, dim=c(iter,Ks))
  nobs <- nrow(X)
  scale_smp <- array(NA, dim=c(iter, nobs))
  lp_smp <- array(NA, dim=c(iter, nobs))
  loglike_smp <- rep(NA, iter)
  logpost_smp <- rep(NA, iter)
  
  # excluding zero-variance variables from sampling process
  idx <- c(1, setdiff(1:K, which(apply(X, 2, function(x) sd(x)<sd.thresh))))
  idxs <- c(1, setdiff(1:Ks, which(apply(Xs, 2, function(x) sd(x)<sd.thresh))))
  if (ordweib) {
    idxs <- 1
  }
  N <- length(idx)
  Ns <- length(idxs)
  
  if (is.null(init)) {
    beta <- rep(0, K)
    betas <- rep(0, Ks)
  } else {
    beta <- init$beta; beta[-idx] <- 0
    betas <- init$betas; betas[-idxs] <- 0
  }
  
  for (n in 1:iter) {
    #bsgw.logpost <- function(coeff, X, Xs, t, s, lambda, lambdas, ...)
    coeffs <- MfU.Sample(c(beta[idx], betas[idxs]), f=bsgw.logpost, uni.sampler="slice"
      , X=X[,idx,drop=FALSE], Xs=Xs[,idxs,drop=FALSE], t=t, s=s, lambda=lambda, lambdas=lambdas
      , ordweib=ordweib, alpha.min=alpha.min, alpha.max=alpha.max
      , control=MfU.Control(n=N+Ns, slice.w=w, slice.m=mbeta
        , slice.lower=c(c(-Inf,rep(-beta.max,N-1)),c(-Inf,rep(-betas.max,Ns-1)))
        , slice.upper=c(c(+Inf,rep(+beta.max,N-1)),c(+Inf,rep(+betas.max,Ns-1)))
        )
      )
    #coeffs <- bsgw.multislice.from.unislice(c(beta[idx], betas[idxs]), bsgw.logpost, X[,idx,drop=F], Xs[,idxs,drop=F]
    #                                             , t, s, lambda, lambdas
    #                                             , ordweib=ordweib, alpha.min=alpha.min, alpha.max=alpha.max, w=w, m=mbeta
    #                                        , lower=c(c(-Inf,rep(-beta.max,N-1)),c(-Inf,rep(-betas.max,Ns-1)))
    #                                        , upper=c(c(+Inf,rep(+beta.max,N-1)),c(+Inf,rep(+betas.max,Ns-1)))
    #                                        )
    beta[idx] <- coeffs[1:N]
    betas[idxs] <- coeffs[N+1:Ns]
    beta_smp[n,] <- beta
    betas_smp[n,] <- betas
    
    ret <- bsgw.eval(coeffs, X, Xs, t, ordweib)
    scale_smp[n,] <- ret$scale
    lp_smp[n,] <- ret$lp
    loglike_smp[n] <- bsgw.loglike(coeffs, X[,idx,drop=F], Xs[,idxs,drop=F], t, s, ordweib, alpha.min, alpha.max)
    logpost_smp[n] <- bsgw.logpost(coeffs, X[,idx,drop=F], Xs[,idxs,drop=F], t, s, lambda, lambdas, ordweib, alpha.min, alpha.max)
    
    if (n%%nskip==0 && print.level>=2) cat("finished sample", n, "of", iter, "\n")
    flush.console()
  }
  return (list(beta=beta_smp, betas=betas_smp, scale=scale_smp, lp=lp_smp, loglike=loglike_smp, logpost=logpost_smp, idx=idx, idxs=idxs))
}



