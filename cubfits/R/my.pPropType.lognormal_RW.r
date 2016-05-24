### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw Gibbs Sampler given current status for measure error (sigmaW) and
### drew inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropType.lognormal_RW <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  log.phi.Curr <- log(phi.Curr)
  p.DrawScale <- .cubfitsEnv$all.DrawScale$p[1]
  p.DrawScale.prev <- .cubfitsEnv$all.DrawScale$p.prev[1]

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


  ### Propose sigma.Phi.Curr.
  proplist <- my.proposesigmaPhi.RW_Norm(
                sigma.Phi.Curr,
                sigma.Phi.DrawScale = p.DrawScale,
                sigma.Phi.DrawScale.prev = p.DrawScale.prev)

  ### M-H step.
  list.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.drawRestrictHP(proplist, list.Curr, phi.Curr)

  ### Update prior's acceptance and adaptive.
  my.update.acceptance("p", ret$accept)
  my.update.adaptive("p", ret$accept)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(sigmaW.Curr, ret$nu.Phi, ret$sigma.Phi)
  ret
} # my.pPropType.lognormal_RW().
