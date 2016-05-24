## 1: make
make.lambda <- function(tree, states, states.sd=0, control=list()) {
  control <- check.control.continuous(control)
  cache <- make.cache.lambda(tree, states, states.sd, control)

  if (control$method == "vcv") {
    all.branches <- make.all.branches.rescale.vcv(cache, control)
    rootfunc <- rootfunc.bm.vcv
  } else if (control$method == "pruning") {
    all.branches <- make.all.branches.lambda.pruning(cache, control)
    rootfunc <- rootfunc.bm.pruning
  } else if (control$method == "contrasts") {
    all.branches <- make.all.branches.rescale.contrasts(cache, control)
    rootfunc <- rootfunc.bm.contrasts
  } else {
    stop("Unknown method", control$method)
  }

  ll <- function(pars, root=ROOT.MAX, root.x=NULL,
                 intermediates=FALSE) {
    check.pars.lambda(pars)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, root, root.x, intermediates)
  }
  class(ll) <- c("lambda", "dtlik", "function")
  ll
}

## 2: info
make.info.lambda <- function(phy) {
  list(name="lambda",
       name.pretty="Pagel's lambda",
       ## Parameters:
       np=2L,
       argnames=default.argnames.lambda(),
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
         "Pagel (1999)"))
}
default.argnames.lambda <- function()
  c("s2", "lambda")

## 3: make.cache
make.cache.lambda <- function(tree, states, states.sd, control) {
  cache <- make.cache.bm(tree, states, states.sd, control)
  cache$info <- make.info.lambda(tree)
  cache
}

###########################################################################
## Additional functions

## Checking
check.pars.lambda <- function(pars) {
  if (length(pars) != 2)
    stop("Incorrect parameter length")
  check.nonnegative(pars)
  ## Technically this is over strict, but real value is complicated.
  if (pars[[2]] > 1)
    stop("lambda must be in [0,1]")
  TRUE
}

make.all.branches.lambda.pruning <- function(cache, control) {
  ## NOTE: This is a hack, but allow here for the extra parameters
  cache$info$np <- 4L

  pars.extra <- c(max(cache$depth), cache$n.tip)
  
  if (control$backend == "R") {
    all.branches <- function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache,
                          initial.conditions.bm.pruning,
                          branches.lambda, preset)
  } else {
    all.branches <- make.all.branches.continuous(cache, control)
  }
  function(pars, ...)
    all.branches(c(pars, pars.extra), ...)
}

## I'm not very happy with the extra information coming in here as
## parameters, but it's the only way that I can make this work easily.
## For the R backend alone it's easy enough to get the extra
## parameters found through a closure (so, turn this into
## make.branches.lambda as we only need two numbers that are constant
## over all parameters.  So unless the "continuous" approach is
## modified to take extra parameters, this is going to be ugly always.
branches.lambda <- function(y, len, pars, t0, idx) {
  m <- y[[1]]
  v <- y[[2]]
  z <- y[[3]]

  sigma2 <- pars[[1]]
  lambda <- pars[[2]]
  tr     <- pars[[3]] # same as eb
  n.tip  <- pars[[4]] # extra, extra, read all about it.

  if (idx > n.tip)
    len.scaled <- len * lambda
  else
    len.scaled <- len * lambda + (1 - lambda) * (tr - t0)

  list(z, c(m, v + sigma2 * len.scaled, 0))
}
