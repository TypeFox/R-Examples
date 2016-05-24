# Contributed by Jonathan Rolland <Jonathan Rolland
# <jonathan.rolland@yahoo.fr>
# based on model-bisse-t.R, and modified by RGF
make.geosse.t <- function(tree, states, functions, sampling.f=NULL,
                          strict=TRUE, control=list(),
                          truncate=FALSE, spline.data=NULL) {
  cache <- make.cache.geosse.t(tree, states,  functions, sampling.f,
                               strict, truncate, spline.data)

  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.geosse)

   rootfunc <- make.rootfunc.t(cache, rootfunc.geosse)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("geosse.t", "geosse", "dtlik", "function")
  ll
}

make.cache.geosse.t <- function(tree, states, functions,
                                sampling.f, strict,
                                truncate, spline.data) {
  cache <- make.cache.geosse(tree, states,sampling.f, strict)
  update.cache.t(cache, functions,
                 nonnegative=TRUE, truncate=truncate, with.q=FALSE,
                 spline.data=spline.data)
}
