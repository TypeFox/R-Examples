bayesmixsurv.eval <- function(alpha1, beta1, alpha2, beta2, X1, X2, gamma, Xg, t) {
  lp1 <- X1%*%beta1
  lp2 <- X2%*%beta2
  if (!is.null(Xg)) {
    lp1 <- lp1 + Xg%*%gamma
    lp2 <- lp2 + Xg%*%gamma
  }
  expterm1 <- exp(lp1)
  expterm2 <- exp(lp2)
  H1 <- t^(alpha1) * expterm1
  H2 <- t^(alpha2) * expterm2
  H <- H1+H2
  h1 <- alpha1 * t^(alpha1-1) * expterm1
  h2 <- alpha2 * t^(alpha2-1) * expterm2
  h <- h1+h2
  S <- exp(-H)
  return (list(h1=h1, h2=h2, h=h, H1=H1, H2=H2, H=H, S=S, lp1=lp1, lp2=lp2))
}

condprob.data <- function(alpha1, beta1, alpha2, beta2, X1, X2, gamma, Xg, t, s) {
  ret <- bayesmixsurv.eval(alpha1, beta1, alpha2, beta2, X1, X2, gamma, Xg, t)
  return (sum((s*log(ret$h)-ret$H)))
}

condprob.beta <- function(lambda, beta) {
  return (-lambda*sum(abs(beta)))
}

condprob.gamma <- function(mu.gamma, sigma.gamma, gamma) {
  return (-length(gamma*log(sigma.gamma)-0.5*sum((gamma-mu.gamma)^2)/sigma.gamma^2))
}

condprob.alpha <- function(alpha) {
  log(dgamma(alpha, shape=1.0, rate=1e-3))
}

## posterior distributions for stochastic nodes ##
condpost.beta <- function(beta.this, beta.that, alpha.this, alpha.that, gamma, X.this, X.that, t, s, lambda.this, Xg) {
  return (condprob.data(alpha.this, beta.this, alpha.that, beta.that, X.this, X.that, gamma, Xg, t, s) + condprob.beta(lambda.this, beta.this))
}
condpost.alpha <- function(alpha.this, beta.this, beta.that, alpha.that, gamma, X.this, X.that, t, s, Xg) {
  return (condprob.data(alpha.this, beta.this, alpha.that, beta.that, X.this, X.that, gamma, Xg, t, s) + condprob.alpha(alpha.this))
}
condpost.gamma <- function(gamma, beta1, beta2, alpha1, alpha2, mu.gamma, sigma.gamma, X1, X2, t, s, Xg) {
  return (condprob.gamma(gamma, mu.gamma, sigma.gamma) + condprob.data(alpha1, beta1, alpha2, beta2, X1, X2, gamma, Xg, t, s))
}

# MCMC routine
bayesmixsurv.mcmc <- function(X1, X2, t, s, Xg, lambda1=1.0, lambda2=lambda1
                              , iter=1000, single=FALSE, alpha2.fixed=NULL, alpha.boundary=1.0
                              , sd.thresh=1e-04, print.level=2, nskip=round(iter/10)) {
  # checking of first column of X1 and X2 being intercept must happen in caller (bayesmixsurv)
  if (single) {
    alpha2.fixed <- 1.0
  }
  
  if (is.null(alpha2.fixed)) {
    alpha1.lower <- 0.001
    alpha1.upper <- alpha.boundary-0.001
    alpha2.lower <- alpha.boundary+0.001
    alpha2.upper <- 10.0 # turn into parameter, perhaps hidden from user
  } else {
    alpha1.lower <- 0.001
    alpha1.upper <- 10.0
  }
  
  K1 <- ncol(X1)
  K2 <- ncol(X2)
  if (!is.null(Xg)) {
    Ncenter <- ncol(Xg)
  } else {
    Ncenter <- 2 # check that 1 works too and switch to 1
  }
  gamma.smp <- array(NA, dim=c(iter,Ncenter))
  beta1.smp <- array(NA, dim=c(iter,K1))
  beta2.smp <- array(NA, dim=c(iter,K2))
  alpha1.smp <- rep(NA, iter)
  alpha2.smp <- rep(NA, iter)
  #mu.gamma.smp <- rep(NA, iter) # leftover code, but mu.gamma is really zero since it represents difference between centers
  sigma.gamma.smp <- rep(NA, iter)
  loglike.smp <- rep(NA, iter)
  
  # excluding zero-variance variables from sampling process
  idx1 <- c(1, setdiff(1:K1, which(apply(X1, 2, function(x) sd(x)<sd.thresh)))) # assumes intercept is at 1; TODO: we need to ensure that beta1[1] corresponds to intercept
  idx2 <- c(1, setdiff(1:K2, which(apply(X2, 2, function(x) sd(x)<sd.thresh)))) # assumes intercept is at 1; TODO: we need to ensure that beta1[1] corresponds to intercept
  
  # do we need better initialization, and allow for inits to be passed in?
  beta1 <- rep(0, K1)
  beta2 <- rep(0, K2)
  if (is.null(alpha2.fixed)) {
    alpha1 <- 0.5
    alpha2 <- 2.0
  } else {
    alpha1 <- 1.0 # can we do better? be more intelligent?
    alpha2 <- alpha2.fixed
  }
  gamma <- rep(0, Ncenter)
  mu.gamma <- 0.0
  sigma.gamma <- 1.0
  
  for (n in 1:iter) {
    # beta's
    beta1[idx1] <- bayesmixsurv.multislice.from.unislice(beta1[idx1], condpost.beta, beta2, alpha1, alpha2, gamma, X1[,idx1,drop=F], X2, t, s, lambda1, Xg
                                            , w=1.0, m=100
                                            , lower=rep(-Inf,length(idx1)), upper=rep(+Inf,length(idx1))) # add maximum absolute value on beta1
    if (!single) {
      beta2[idx2] <- bayesmixsurv.multislice.from.unislice(beta2[idx2], condpost.beta, beta1, alpha2, alpha1, gamma, X2[,idx2,drop=F], X1, t, s, lambda2, Xg
                                              , w=1.0, m=100
                                              , lower=rep(-Inf,length(idx2)), upper=rep(+Inf,length(idx2)))
    } else {
      beta2[1] <- -1e+12 # forcing e^(xbeta2) to be 0; TODO: we need to ensure that beta1[1] corresponds to intercept
    }

    # alpha's
    alpha1 <- bayesmixsurv.uni.slice(alpha1, condpost.alpha, beta1, beta2, alpha2, gamma, X1, X2, t, s, Xg
                        , w=1.0, m=10, lower=alpha1.lower, upper=alpha1.upper)
    if (is.null(alpha2.fixed)) {
      alpha2 <- bayesmixsurv.uni.slice(alpha2, condpost.alpha, beta2, beta1, alpha1, gamma, X2, X1, t, s, Xg
                          , w=1.0, m=10, lower=alpha2.lower, upper=alpha2.upper)
    }

    # gamma
    if (!is.null(Xg)) {
      gamma <- bayesmixsurv.multislice.from.unislice(gamma, condpost.gamma, beta1, beta2, alpha1, alpha2, mu.gamma, sigma.gamma, X1, X2, t, s, Xg
                                        , w=1.0, m=10, lower=-Inf, upper=+Inf)
      sigma.gamma <- sqrt(1/rgamma(1, shape=Ncenter/2, scale=0.5*sum((mu.gamma-gamma)^2))) # using conjugacy, what about prior?! did i assume flat prior and just ignore it?
    }
    
    beta1.smp[n,] <- beta1
    beta2.smp[n,] <- beta2
    alpha1.smp[n] <- alpha1
    alpha2.smp[n] <- alpha2
    gamma.smp[n,] <- gamma
    #mu.gamma.smp[n] <- mu.gamma
    sigma.gamma.smp[n] <- sigma.gamma
    
    loglike.smp[n] <- condprob.data(alpha1, beta1, alpha2, beta2, X1, X2, gamma, Xg, t, s)
    
    if (n%%nskip==0 && print.level>=2) cat("finished sample", n, "of", iter, "\n")
    flush.console() # for windows
  }
  return (list(beta1=beta1.smp, beta2=beta2.smp, alpha1=alpha1.smp, alpha2=alpha2.smp, loglike=loglike.smp
               , gamma=gamma.smp
               , sigma.gamma=sigma.gamma.smp#, xbeta1=xbeta1, xbeta2=xbeta2
               , idx1=idx1, idx2=idx2))
}

