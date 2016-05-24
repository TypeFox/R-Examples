### Gibbs sampler and M-H for step 2, mainly for hyperparameters.

### Draw prior for log mixture normal distribtuion.
my.pPropType.logmixture <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  paramlog.Curr <- p.Curr[-1]
  log.phi.Curr <- log(phi.Curr)
  p.DrawScale <- .cubfitsEnv$all.DrawScale$p
  p.DrawScale.prev <- .cubfitsEnv$all.DrawScale$p.prev

  ### Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  ###                            (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
  
  if(.CF.CONF$estimate.Phi.noise)
  {
    sigmaW.Curr <-
      sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                      rate = sum((log.phi.Obs - log.phi.Curr)^2) / 2))
  }else{
    sigmaW.Curr <- p.Curr[1]
  }


  ### Propose new mixture parameters.
  proplist <- my.propose.paramlog(paramlog.Curr, log.phi.Curr, hp.param)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", 1)
  my.update.adaptive("p", 1)

  ### Return.
  ret <- c(sigmaW.Curr, proplist$paramlog)
  ret
} # my.pPropType.logmixture().


### Ref: Bayesian Modelling and Inference on Mixtures of Distributions,
###      Jean-Michel Marin, Kerrie Mengersen and Christian P. Robert,
###      page 31.
### Ref: https://www.ceremade.dauphine.fr/~xian/BCS/

### Propose new param in log scale for log mixture normal distribution.
my.propose.paramlog <- function(paramlog.Curr, log.phi.Curr, hp.param){
  ### Dispatch.
  K <- length(hp.param$hp.Init) / 3
  gamma.k <- hp.param$hp.Init[1:K] + 1    # + 1 to avoid degeneration.
  lambda.k <- hp.param$hp.Init[(1:K) + K]
  sigma.k <- hp.param$hp.Init[(1:K) + K + K]
  sigma.Curr <- paramlog.Curr[(1:K) + K + K]

  ### Preset hyperparameters.
  ### sigma.k^2 ~ IG(alpha.k, 1 / beta.k) for k = 1, ..., K
  ### mu.k | sigma.k^2 ~ N(lambda.k, sigma.k^2 / tau.k) for k = 1, ..., K
  ### eta ~ D(gamma.1, ..., gamma.K)
  alpha.k <- rep(1.0, K)
  beta.k <- sigma.k
  tau.k <- rep(1.0, K)

  ### Draw missing indicators.
  ### Z \propto eta.k f(x_i | theta.k)
  Z.prop <- proplmixnorm(log.phi.Curr, paramlog.Curr)
  Z <- apply(Z.prop, 1, function(x){ sample.int(K, size = 1, prob = x) })

  ### Draw eta.k from conjugates.
  ### eta ~ D(r0.1 + n.1^{(t)}, ..., r0.K + n.K^{(t)})
  ### n.k = sum_{i=1}^N I(z_i == k)
  n.k <- as.vector(tabulate(Z, nbins = K))
  eta.New <- rdirichlet(1, n.k + gamma.k)

  ### Draw mu.k from conjugates.
  ### mu.k ~ N((lambda.k * tau.k + s.m.k^{(t)}) / (tau.k + n.k^{(t)}),
  ###          sigma.k^2^{(t - 1)} / (tau.k + n.k^{(t)})
  ### s.m.k^{(t)} = sum_{i=1}^N I(z_i == k) x_i
  mu.New <- rep(0.0, K)
  for(i.k in 1:K){
    s.m.k <- sum(log.phi.Curr[Z == i.k])
    mean.New <- (lambda.k[i.k] * tau.k[i.k] + s.m.k) / (tau.k[i.k] + n.k[i.k])
    sd.New <- sigma.Curr[i.k] / sqrt(tau.k[i.k] + n.k[i.k])
    mu.New[i.k] <- rnorm(1, mean.New, sd = sd.New)
  }

  ### Draw sigma.k from conjugates.
  ### sigma.k ~ IG(alpha.k + (n.k + 1) / 2,
  ###              (beta.k + tau.k * (mu.k^{(t)} - lambda.k)^2 / 2 +
  ###                        s.v.k^{(t)} / 2))
  ### s.v.k^{(t)} = sum_{i=1}^N I(z_i == k) (x_i - mu.k^{(t)})^2
  sigma.New <- rep(0.0, K)
  for(i.k in 1:K){
    s.v.k <- sum((log.phi.Curr[Z == i.k] - mu.New[i.k])^2)
    alpha.New <- alpha.k[i.k] + (n.k[i.k] + 1) * 0.5
    beta.New <- beta.k[i.k] +
                0.5 * tau.k[i.k] * (mu.New[i.k] - lambda.k[i.k])^2 +
                0.5 * s.v.k
    sigma.New[i.k] <- 1 / rgamma(1, shape = alpha.New, rate = beta.New)
  }
  sigma.New <- sqrt(sigma.New)

  ### Return.
  ret <- list(paramlog = as.numeric(c(eta.New, mu.New, sigma.New)))
  ret
} # End of my.propose.paramlog().
