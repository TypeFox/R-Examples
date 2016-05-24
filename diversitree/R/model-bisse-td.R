make.bisse.td <- function(tree, states, n.epoch, unresolved=NULL,
                              sampling.f=NULL, nt.extra=10,
                              strict=TRUE, control=list()) {
  cache <- make.cache.bisse.td(tree, states, n.epoch, unresolved,
                               sampling.f, nt.extra, strict)
  all.branches <- make.all.branches.td.dtlik(cache, control,
                                             initial.conditions.bisse)
  rootfunc <- make.rootfunc.td(cache, rootfunc.musse)
  f.pars <- make.pars.bisse.td(n.epoch)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.td", "bisse", "dtlik", "function")
  ll
}

## 3: argnames / argnames<-
make.cache.bisse.td <- function(tree, states, n.epoch, unresolved,
                                sampling.f, nt.extra, strict) {
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying BiSSE with unresolved clades")
  cache$info <- update.info.td(cache$info, n.epoch)
  cache
}

make.pars.bisse.td <- function(n.epoch) {
  npar <- (n.epoch - 1) + (6 * n.epoch)
  i.t <- seq_len(n.epoch - 1)
  i.p <- n.epoch:npar
  function(pars) {
    check.pars.nonnegative(pars, npar)
    cbind(c(pars[i.t], Inf), matrix(pars[i.p], n.epoch, 6, TRUE))
  }
}
