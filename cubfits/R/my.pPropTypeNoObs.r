### Gibbs sampler for step 2, mainly for hyperparameters.

### Get the specific function according to the options.
get.my.pPropTypeNoObs <- function(type){
  if(!any(type[1] %in% .CF.CT$type.p)){
    stop("type is not found.")
  }
  ret <- eval(parse(text = paste("my.pPropTypeNoObs.", type[1], sep = "")))
  assign("my.pPropTypeNoObs", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.pPropTypeNoObs().


### Draw Gibbs Sampler given current status for lognormal prior around fixed
### mean of log expression.
my.pPropTypeNoObs.lognormal_fix <- function(n.G, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[1]
  sigma.Phi.Curr <- p.Curr[2]
  log.phi.Obs.mean <- hp.param$log.phi.Obs.mean
  log.phi.Curr <- log(phi.Curr)

  ### nu.Phi.Curr is unused in this case.
  ### sigma.Phi.Curr is unused in this case.

  ### The next is assuming non-informative priors for mu.Phi and sigma.Phi.sq.
  ### i.e. p(mu.Phi, sigma.Phi.sq) \propto 1/sigma.Phi.sq.
  ### See Gelman et al. (2003), p.75 for details.

  ### Draw \sigma^{2*}_phi from IG((n_G - 1) / 2, (n_G - 1) S^{2(t)}_{phi} / 2)
  ###   sigma.Phi.Curr <- sqrt(1 / rgamma(1, shape = (n_G - 1) / 2,
  ###                                        rate = sum((log.phi.Curr -
  ###                                                    nu.Phi.Curr)^2) / 2))
  ### Draw \mu^*_phi from N(sum_g log(phi^{(t)}_g) / n_G, \sigma^{2*}_phi / n_G)
  # nu.Phi.Curr <- mean(log.phi.Curr) + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  ### Since there is no phi.Obs, log.phi.Curr can draft potentially to anywhere.
  ### This is to avoid draft of Phi where log.phi.Obs.mean is fixed for all
  ### iterations, but we still allow some uncertainty around the mean
  ### value for all Phi's.
  sigma.Phi.Curr <-
    sqrt(1 / rgamma(1, shape = (n.G - 1) / 2,
                       rate = sum((log.phi.Curr - log.phi.Obs.mean)^2) / 2))
  nu.Phi.Curr <- log.phi.Obs.mean + rnorm(1) * sigma.Phi.Curr / sqrt(n.G)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", 1)
  my.update.adaptive("p", 1)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(nu.Phi.Curr, sigma.Phi.Curr)
  ret
} # End of my.pPropTypeNoObs.lognormal_fix().

