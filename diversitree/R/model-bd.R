## The B-D model is different to the others in that I am not using
## most of the infrastructure - instead the entire calculations are
## done at once.

## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches
make.bd <- function(tree, sampling.f=NULL, unresolved=NULL,
                    times=NULL, control=list()) {
  control <- check.control.bd(control, times)
  cache <- make.cache.bd(tree, sampling.f, unresolved, times, control)
  const <- cache$const

  if ( control$method == "nee" ) {
    all.branches <- make.all.branches.bd.nee(cache, control)
    rootfunc <- rootfunc.bd.nee
  } else {
    all.branches <- make.all.branches.dtlik(cache, control,
                                            initial.conditions.bd.ode)
    rootfunc <- rootfunc.bd.ode
  }

  ll <- function(pars, condition.surv=TRUE, intermediates=FALSE) {
    check.pars.nonnegative(pars, 2)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, intermediates, const)
  }
  class(ll) <- c("bd", "dtlik", "function")
  ll
}

## Yule is somewhat weird, as we allow likelihood calculations, but
## cheat on the ML search and go straight for the ML point.  Well, we
## used to, but that is currently broken.
make.yule <- function(tree, sampling.f=NULL, unresolved=NULL,
                      times=NULL, control=list()) {
  control <- check.control.bd(control)
  ll.bd <- make.bd(tree, sampling.f, unresolved, times, control)
  ll <- function(pars, condition.surv=TRUE)
    ll.bd(c(pars, 0), condition.surv)
  set.info(ll, make.info.yule(tree))
  class(ll) <- c("yule", "bd", "dtlik", "function")
  ll
}

## 2: info
## TODO: This should deal with the case where 'times' but not 'phy' is
## given?
make.info.bd <- function(phy) {
  list(name="bd",
       name.pretty="Constant rate birth-death",
       ## Parameters:
       np=2L,
       argnames=default.argnames.bd(),
       ## Variables:
       ny=2L,
       k=1L,
       idx.e=1L,
       idx.d=2L,
       ## R version of the derivatives function
       derivs=derivs.bd,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="nlm",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Nee et al. (1994) doi:10.1098/rstb.1994.0068"))
}
default.argnames.bd <- function()
  c("lambda", "mu")

make.info.yule <- function(phy) {
  info <- make.info.bd(phy)
  info$name.yule <- "yule"
  info$name.pretty <- "Yule (pure birth)"
  info$argnames <- default.argnames.yule()
  info
}
default.argnames.yule <- function()
  "lambda"

make.cache.bd <- function(tree, sampling.f, unresolved, times,
                          control) {
  if ( control$method == "nee" )
    cache <- make.cache.bd.nee(tree, sampling.f, unresolved, times)
  else
    cache <- make.cache.bd.ode(tree, sampling.f, unresolved)
  cache$info <- make.info.bd(tree)
  cache
}

######################################################################
## Extra functions
starting.point.bd <- function(tree, yule=FALSE) {
  p.yule <- c(stats::coef(suppressWarnings(find.mle(make.yule(tree), .1))), 0)
  names(p.yule) <- default.argnames.bd()
  if (yule)
    p.yule
  else
    suppressWarnings(stats::coef(find.mle(make.bd(tree), p.yule)))
}

## Checking
check.control.bd <- function(control, times) {
  control <- modifyList(list(method="nee"), control)
  if ( !(control$method %in% c("nee", "ode")) )
    stop(sprintf("Unknown method '%s'", control$method))
  control
}

check.unresolved.bd <- function(tree, unresolved) {
  if ( is.null(tree) )
    stop("Cannot just specify times when using unresolved clades")
  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- make.unresolved.bd(tree$clades)
  } else if ( !is.null(unresolved) && length(unresolved) > 0 ) {
    if ( is.null(names(unresolved)) || !is.numeric(unresolved) )
      stop("'unresolved' must be a named numeric vector")
    if ( !(all(names(unresolved) %in% tree$tip.label)) )
      stop("Unknown species in 'unresolved'")
    if ( any(unresolved < 1) )
      stop("All unresolved entries must be > 0")

    if ( any(unresolved == 1) ) {
      warning("Removing unresolved entries that are one")
      unresolved <- unresolved[unresolved != 1]
    }
  }

  if ( length(unresolved) == 0 )
    unresolved <- NULL
  else {
    i <- match(names(unresolved), tree$tip.label)
    list(n=unresolved,
         t=tree$edge.length[match(i, tree$edge[,2])])
  }
}
