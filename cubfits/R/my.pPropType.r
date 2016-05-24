### Gibbs sampler for step 2, mainly for hyperparameters.

### Get the specific function according to the options.
get.my.pPropType <- function(type){
  if(!any(type[1] %in% .CF.CT$type.p)){
    stop("type is not found.")
  }
  ret <- eval(parse(text = paste("my.pPropType.", type[1], sep = "")))
  assign("my.pPropType", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.pPropType().


### Draw Gibbs Sampler given current status for lognormal prior around fixed
### mean of log expression.
my.pPropType.lognormal_fix <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  log.phi.Obs.mean <- hp.param$log.phi.Obs.mean
  log.phi.Curr <- log(phi.Curr)

  ### sigma.Phi.Curr is unused in this case.

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

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ### As my.pPropTypeNoObs.lognormal_fix(), we can set a fixed mean to
  ### log.phi.Obs if phi.Obs has been normalized to mean=1 prior.
  ### Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.phi.Curr - log.phi.Obs.mean)^2) / 2))
  ### Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  nu.Phi.Curr <- log.phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", 1)
  my.update.adaptive("p", 1)

  ### Return.
  ret <- c(sigmaW.Curr, nu.Phi.Curr, sigma.Phi.Curr)
  ret
} # my.pPropType.lognormal_fix().

