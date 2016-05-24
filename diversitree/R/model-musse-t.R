make.musse.t <- function(tree, states, k, functions, sampling.f=NULL,
                         strict=TRUE, control=list(),
                         truncate=FALSE, spline.data=NULL) {
  cache <- make.cache.musse.t(tree, states, k, functions, sampling.f,
                              strict, truncate, spline.data)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.musse)
  rootfunc <- make.rootfunc.t(cache, rootfunc.musse)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("musse.t", "musse", "dtlik.t", "dtlik", "function")
  ll
}

make.cache.musse.t <- function(tree, states, k, functions,
                               sampling.f, strict,
                               truncate, spline.data) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  update.cache.t(cache, functions,
                 nonnegative=TRUE, truncate=truncate, with.q=TRUE,
                 spline.data=spline.data) 
}
