### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw Gibbs Sampler given current status for measure error (sigmaW) and
### drew inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression. Also, draw bias.Phi.
my.pPropType.lognormal_bias <- function(n.G, log.phi.Obs, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[2]
  sigma.Phi.Curr <- p.Curr[3]
  bias.Phi.Curr <- p.Curr[4]
  log.phi.Curr <- log(phi.Curr)
  p.DrawScale <- .cubfitsEnv$all.DrawScale$p
  p.DrawScale.prev <- .cubfitsEnv$all.DrawScale$p.prev

  ### Draw \sigma^{2*}_W from IG((n_G - 1) / 2,
  ###                         (n_G - 1) S^{2(t)}_{phi_{obs}} / 2)
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
                sigma.Phi.DrawScale = p.DrawScale[1],
                sigma.Phi.DrawScale.prev = p.DrawScale.prev[1])

  ### M-H step for hyperparameters.
  list.Curr <- list(nu.Phi = nu.Phi.Curr, sigma.Phi = sigma.Phi.Curr)
  ret <- my.drawRestrictHP(proplist, list.Curr, phi.Curr)

  ### Propose bias.Phi.
  proplist <- my.proposebiasPhi.RW_Norm(
                bias.Phi.Curr,
                bias.Phi.DrawScale = p.DrawScale[2],
                bias.Phi.DrawScale.prev = p.DrawScale.prev[2])

  ### M-H step for hyperparameters of bias.
  list.Curr <- list(bias.Phi = bias.Phi.Curr,
                    nu.Phi = ret$nu.Phi, sigma.Phi = ret$sigma.Phi)
  ### ret$nu.Phi and ret$sigma.Phi are in common of bias.Phi.Curr and
  ### proplist$bias.Phi, so there is no need to make a copy in proplist.
  ### However, both ret$nu.Phi and ret$sigma.Phi should be returned at the end.
  ret.bias <- my.drawbiasPhi(proplist, list.Curr, log.phi.Obs, log.phi.Curr,
                             sigmaW.Curr)

  ### Update prior's acceptance and adaptive.
  accept <- c(ret$accept, ret.bias$accept)
  my.update.acceptance("p", accept)
  my.update.adaptive("p", accept)

  ### Only nu.Phi and sigma.Phi are used.
  ret <- c(sigmaW.Curr, ret$nu.Phi, ret$sigma.Phi, ret.bias$bias.Phi)
  ret
} # my.pPropType.lognormal_bias().
