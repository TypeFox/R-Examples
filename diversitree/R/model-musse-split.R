## Split models should provide
##   1. make
##   2. make.cache
## Common other functions include:
##   branches.aux

## 1: make
make.musse.split <- function(tree, states, k, nodes, split.t,
                             sampling.f=NULL, strict=TRUE,
                             control=list()) {
  cache <- make.cache.musse.split(tree, states, k, nodes, split.t,
                                  sampling.f, strict)
  n.part <- cache$n.part

  all.branches <- make.all.branches.split.dtlik(cache, control,
                                                initial.conditions.musse)
  rootfunc <- make.rootfunc.split(cache, rootfunc.musse)
  f.pars <- make.pars.musse.split(k, n.part)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }
  
  class(ll) <- c("musse.split", "musse", "dtlik", "function")
  ll
}

## Make requires the usual functions:
## 2: make.cache (initial.tip, root)
make.cache.musse.split <- function(tree, states, k, nodes, split.t,
                                   sampling.f, strict) {
  cache <- make.cache.musse(tree, states, k, NULL, strict)
  make.cache.split.xxsse(tree, cache, nodes, split.t, sampling.f)
}

make.branches.aux.musse <- function(cache, control)
  make.branches.aux.dtlik(cache, control)

make.pars.musse.split <- function(k, n.part) {
  f.pars <- make.pars.musse(k)
  function(pars) {
    pars <- check.pars.multipart(pars, n.part, k*(k+1))
    for ( i in seq_len(n.part) )
      pars[[i]] <- f.pars(pars[[i]])
    pars
  }
}
