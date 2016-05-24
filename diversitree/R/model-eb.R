## My simple-minded EB calculator.

## This is a bit of a trick, because we need to augment everything
## with the depth.

## 1: make
make.eb <- function(tree, states, states.sd=0, control=list()) {
  control <- check.control.continuous(control)
  cache <- make.cache.eb(tree, states, states.sd, control)

  if (control$method == "vcv") {
    all.branches <- make.all.branches.rescale.vcv(cache, control)
    rootfunc <- rootfunc.bm.vcv
  } else {
    all.branches <- make.all.branches.eb.pruning(cache, control)
    rootfunc <- rootfunc.bm.pruning
  }

  ll <- function(pars, root=ROOT.MAX, root.x=NULL,
                 intermediates=FALSE) {
    check.pars.eb(pars)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, root, root.x, intermediates)
  }
  class(ll) <- c("eb", "dtlik", "function")
  ll
}

## 2: info
make.info.eb <- function(phy) {
  list(name="eb",
       name.pretty="Early Burst (AC/DC)",
       ## Parameters:
       np=2L,
       argnames=default.argnames.eb(),
       ## Variables:
       ny=3L,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=FALSE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Blomberg et al. 2003. Evolution 57: 717"))
}
default.argnames.eb <- function()
  c("s2", "a")

## 3: make.cache
make.cache.eb <- function(tree, states, states.sd, control) {
  cache <- make.cache.bm(tree, states, states.sd, control)
  cache$info <- make.info.eb(tree)
  cache
}

###########################################################################
## Additional functions

## Checking
check.pars.eb <- function(pars) {
  if (length(pars) != 2)
    stop("Incorrect parameter length")
  check.nonnegative(pars[1])
  check.nonpositive(pars[2])
  TRUE
}

make.all.branches.eb.pruning <- function(cache, control) {
  ## NOTE: This is a hack, but allow here for the extra paramter:
  cache$info$np <- 3L

  pars.extra <- max(cache$depth)

  if (control$backend == "R") {
    all.branches <- function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache,
                          initial.conditions.bm.pruning,
                          branches.eb, preset)
  } else {
    all.branches <- make.all.branches.continuous(cache, control)
  }

  function(pars, ...)
    all.branches(c(pars, pars.extra), ...)
}

## The issue that I have here is that time is computed against the
## root, so we'll need to know that.  t0 is the *depth* of the
## branch *tip* so t0+len is the *depth* of the branch base:
##    tr    t0+len         t0      0
##    |-----|--------------|-------|
##    0     s0             s1      st
## So let, st = tr (time of root, time of tip), then
##    s1 = tr - t0
##    s0 = tr - (t0 + len) = s1 - len
branches.eb <- function(y, len, pars, t0, idx) {
  m <- y[[1]]
  v <- y[[2]]
  z <- y[[3]]

  sigma2 <- pars[[1]]
  a      <- pars[[2]]
  tr     <- pars[[3]]

  if (a != 0) {
    s1 <- tr - t0
    s0 <- s1 - len
    len <- (exp(a * s1) - exp(a * s0))/a
  }

  list(z, c(m, v + sigma2 * len, 0))
}
