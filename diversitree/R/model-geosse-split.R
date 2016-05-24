## GeoSSE model, by Emma Goldberg <eeg@uic.edu>

## Split models should provide
##   1. make
##   2. make.cache
## Common other functions include:
##   branches.aux

## 1: make
make.geosse.split <- function(tree, states, nodes, split.t=Inf,
                              sampling.f=NULL, strict=TRUE,
                              control=list()) {
  cache <- make.cache.geosse.split(tree, states, nodes, split.t,
                                   sampling.f, strict)
  n.part <- cache$n.part

  all.branches <- make.all.branches.split.dtlik(cache, control,
                                                initial.conditions.geosse)
  rootfunc <- make.rootfunc.split(cache, rootfunc.geosse)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.pars.geosse.split(pars, n.part)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
 
  class(ll) <- c("geosse.split", "geosse", "function")
  ll
}

make.geosse.uneven <- function(tree, states, nodes, split.t=Inf,
                             sampling.f=NULL, strict=TRUE,
                             control=list()) {
 cache <- make.cache.geosse.split(tree, states, nodes, split.t,
                                 sampling.f, strict)
 cache$info <- update.info.uneven(cache$info, make.info.geosse(tree))
 n.part <- cache$n.part

 all.branches <- make.all.branches.split.dtlik(cache, control,
                                               initial.conditions.geosse)
 rootfunc <- make.rootfunc.split(cache, rootfunc.geosse)

 ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                root.p=NULL, intermediates=FALSE) {
   pars <- rep(list(check.pars.geosse(pars)), n.part)
   ans <- all.branches(pars, intermediates)
   rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
 }
 class(ll) <- c("geosse.uneven", "geosse", "dtlik", "function")
 ll
}

## Make requires the usual functions:
## 2: make.cache (initial.tip, root)
make.cache.geosse.split <- function(tree, states, nodes, split.t=Inf,
                                    sampling.f, strict) {
  cache <- make.cache.geosse(tree, states, NULL, strict)
  make.cache.split.xxsse(tree, cache, nodes, split.t, sampling.f)
}

make.branches.aux.geosse <- function(cache, control)
  make.branches.aux.dtlik(cache, control)

## when classe.split is working, this can be a special case
check.pars.geosse.split <- function(pars, n.part)
  pars <- check.pars.multipart(check.nonnegative(pars), n.part, 7)
