### Gibbs sampler and random walk for step 2, mainly for hyperparameters.

### Draw inv-/gamma for lognormal priors (nu.Phi, sigma.Phi) around current
### mean of log expression.
my.pPropTypeNoObs.lognormal_RW <- function(n.G, phi.Curr,
    p.Curr, hp.param){
  ### Dispatch.
  nu.Phi.Curr <- p.Curr[1]
  sigma.Phi.Curr <- p.Curr[2]
  p.DrawScale <- .cubfitsEnv$all.DrawScale$p[1]
  p.DrawScale.prev <- .cubfitsEnv$all.DrawScale$p.prev[1]

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
  ret <- c(ret$nu.Phi, ret$sigma.Phi)
  ret
} # my.pPropTypeNoObs.lognormal_RW().
