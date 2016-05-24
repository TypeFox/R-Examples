## 1: make
make.bm <- function(tree, states, states.sd=0, control=list()) {
  control <- check.control.continuous(control)
  cache <- make.cache.bm(tree, states, states.sd, control)

  if (control$method == "vcv") {
    all.branches <- make.all.branches.bm.vcv(cache, control)
    rootfunc <- rootfunc.bm.vcv
  } else if (control$method == "pruning") {
    all.branches <- make.all.branches.bm.pruning(cache, control)
    rootfunc <- rootfunc.bm.pruning
  } else if (control$method == "contrasts") {
    all.branches <- make.all.branches.bm.contrasts(cache, control)
    rootfunc <- rootfunc.bm.contrasts
  } else {
    stop("Unknown method", control$method)
  }

  ll <- function(pars, root=ROOT.MAX, root.x=NULL,
                 intermediates=FALSE) {
    check.pars.nonnegative(pars, 1)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, root, root.x, intermediates)
  }
  class(ll) <- c("bm", "dtlik", "function")
  ll
}

## 2: info
make.info.bm <- function(phy) {
  list(name="bm",
       name.pretty="Brownian motion",
       ## Parameters:
       np=1L,
       argnames=default.argnames.bm(),
       ## Variables:
       ny=3L,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "I really don't know"))
}
default.argnames.bm <- function() "s2"

make.cache.bm <- function(tree, states, states.sd, control) {
  method <- control$method
  tree <- check.tree(tree, ultrametric=FALSE)

  cache <- make.cache(tree)
  if ( is.null(states.sd) || all(states.sd == 0) ) {
    cache$states <- check.states(tree, states, as.integer=FALSE)
    cache$states.sd <- rep(0, cache$n.tip)
  } else {
    tmp <- check.states.quasse(tree, states, states.sd)
    cache$states    <- tmp$states
    cache$states.sd <- tmp$states.sd
  }

  ## There are slightly different things needed for the different
  ## approaches; only compute things where we know that it wil be
  ## useful.
  if (method == "vcv") {
    cache$vcv <- vcv(tree)
  } else if (method == "pruning") {
    cache$y <- initial.tip.bm.pruning(cache)
  } else if (method == "contrasts") {
    ## There is duplication here with make.cache.pgls; perhaps merge
    ## these?  That might help with some of the root treatment
    ## things.
    cache$n <- length(tree$tip.label)
    pics <- pic(cache$states, tree, var.contrasts=TRUE)
    cache$u <- pics[,"contrasts"]
    cache$V <- pics[,"variance"]
    cache$V0 <- pgls.root.var.bm(tree)
    cache$root.x <- pgls.root.mean.bm(tree, cache$states)
  }
  cache$info <- make.info.bm(tree)
  cache
}
