make.bd.t <- function(tree, functions, sampling.f=NULL,
                      unresolved=NULL, control=list(),
                      truncate=FALSE, spline.data=NULL) {
  cache <- make.cache.bd.t(tree, functions, unresolved, sampling.f,
                           truncate, spline.data)
  all.branches <- make.all.branches.t.dtlik(cache, control,
                                            initial.conditions.bd.ode)
  rootfunc <- make.rootfunc.t(cache, rootfunc.bd.ode)
  const <- cache$const
  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, intermediates, const)
  }
  class(ll) <- c("bd.t", "bd", "dtlik.t", "dtlik", "function")
  ll
}

make.cache.bd.t <- function(tree, functions, unresolved, sampling.f,
                            truncate, spline.data) {
  cache <- make.cache.bd.ode(tree, sampling.f, unresolved)
  if ( !is.null(cache$unresolved) )
    stop("Cannot do time-varying birth-death with unresolved clades")
  cache$info$ml.default <- "subplex"
  update.cache.t(cache, functions,
                 nonnegative=TRUE, truncate=truncate, with.q=FALSE,
                 spline.data=spline.data)
}
