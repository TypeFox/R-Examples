## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## A special case of the Mk model where traits are ordered.  This will
## run substantially faster.
##
## I'm making this its own function, rather than a control extension
## to the make.mkn, to parallel make.mkn.multitrait.  Unlike
## mkn.multitrait, this does not actually use the mkn code (much).
make.mkn.meristic <- function(tree, states, k, control=list()) {
  cache <- make.cache.mkn.meristic(tree, states, k, control)
  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions.mkn)
  rootfunc <- rootfunc.mkn

  ll <- function(pars, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) {
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, root, root.p, intermediates)
  }
  class(ll) <- c("mkn.meristic", "mkn", "dtlik", "function")
  ll
}
                              
## 2: info
make.info.mkn.meristic <- function(k, phy) {
  list(name="mkn.meristic",
       name.ode="mkn_meristic",
       name.pretty="Mkn (meristic)",
       ## Parameters:
       np=2L,
       argnames=default.argnames.mkn.meristic(),
       ## Variables:
       ny=as.integer(k),
       k=as.integer(k),
       idx.e=integer(0),
       idx.d=seq_len(k),
       ## R version of derivatives
       derivs=make.derivs.mkn.meristic(k),
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Pagel (1994)",
         "Lewis (2001)"))
}
default.argnames.mkn.meristic <- function(k)
  c("q.down", "q.up")

## 3: cache
make.cache.mkn.meristic <- function(tree, states, k, control) {
  ## Note that in the case of k==2, the root function is slightly
  ## wrong (would need parameters reversed), so I'm just going to
  ## assert that k must be greater than 2.  The more general case is
  ## not handled anyway.
  if ( k <= 2 )
    stop("k must be at least 3 for mkn.meristic")
  
  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=FALSE,
                         strict.vals=1:k)

  cache <- make.cache(tree)
  cache$info <- make.info.mkn.meristic(k, tree)
  cache$states  <- states
  cache$y <- initial.tip.mkn.ode(cache)

  cache
}

## 4: initial conditions [use mkn]
## 5: rootfunc           [use mkn]

mkn.meristic.Q <- function(pars, k) {
  q <- matrix(0, k, k)
  i <- cbind(2:k, 1:(k-1))
  q[i] <- pars[1]
  q[i[,2:1]] <- pars[2]
  diag(q) <- -rowSums(q)
  q
}

make.derivs.mkn.meristic <- function(k) {
  i <- seq(2, k-1)
  dydt <- numeric(k)
  function(t, y, pars) {
    d <- pars[1]
    u <- pars[2]
    ud <- d + u
    dydt[1] <- -u * y[1] + u * y[2]
    dydt[i] <- -ud * y[i] + d * y[i-1] + u * y[i+1]
    dydt[k] <- -d * y[k] + d * y[k-1]
    dydt
  }
}
