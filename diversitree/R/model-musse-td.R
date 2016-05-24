make.musse.td <- function(tree, states, k, n.epoch, sampling.f=NULL, 
                          strict=TRUE, control=list()) {
  cache <- make.cache.musse.td(tree, states, k, n.epoch,
                               sampling.f, strict)
  all.branches <- make.all.branches.td.dtlik(cache, control,
                                             initial.conditions.musse)
  rootfunc <- make.rootfunc.td(cache, rootfunc.musse)
  f.pars <- make.pars.musse.td(n.epoch, k)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("musse.td", "musse", "dtlik", "function")
  ll
}

make.cache.musse.td <- function(tree, states, k, n.epoch, sampling.f,
                                strict) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  cache$info <- update.info.td(cache$info, n.epoch)
  cache
}

make.pars.musse.td <- function(n.epoch, k) {
  np.in <- k * (k + 1)
  np.out <- k * (k + 2) + 1
  npar <- (n.epoch - 1) + (np.in * n.epoch)
  i.t <- seq_len(n.epoch - 1)
  i.p <- n.epoch:npar
  f.pars <- make.pars.musse(k)
  
  function(pars) {
    if ( length(pars) != npar )
      stop(sprintf("Invalid length parameters (expected %d)", npar))

    pars2 <- matrix(NA, n.epoch, np.out)
    pars2[,1] <- c(pars[i.t], Inf)
    tmp <- matrix(pars[i.p], n.epoch, np.in, TRUE)
    for ( i in seq_len(n.epoch) )
      pars2[i,-1] <- f.pars(tmp[i,])
    pars2
  }
}
