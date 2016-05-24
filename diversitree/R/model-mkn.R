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
make.mkn <- function(tree, states, k, strict=TRUE, control=list()) {
  control <- check.control.mkn(control, k)
  cache <- make.cache.mkn(tree, states, k, strict, control)
  all.branches <- make.all.branches.mkn(cache, control)
  rootfunc <- rootfunc.mkn
  f.pars <- make.pars.mkn(k)

  ll <- function(pars, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) {
    qmat <- f.pars(pars)
    ans <- all.branches(qmat, intermediates)
    rootfunc(ans, qmat, root, root.p, intermediates)
  }
  class(ll) <- c("mkn", "dtlik", "function")
  ll
}

make.mk2 <- function(tree, states, strict=TRUE, control=list())
  make.mkn(tree, states + 1, 2, strict, check.control.mk2(control))

## 2: info
make.info.mkn <- function(k, phy) {
  list(name="mkn",
       name.pretty="Mk(n)",
       ## Parameters:
       np=as.integer(k * k),
       argnames=default.argnames.mkn(k),
       ## Variables:
       ny=as.integer(k), # TODO/NEW: only for ode version...
       k=as.integer(k),
       idx.e=integer(0),
       idx.d=seq_len(k),
       ## R version of the derivatives function (only applicable to
       ## ode version, so added there)
       ## derivs=derivs.mkn,
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
default.argnames.mkn <- function(k) {
  base <- ceiling(log10(k + .5))
  fmt <- sprintf("q%%0%dd%%0%dd", base, base)
  sprintf(fmt, rep(1:k, each=k-1),
          unlist(lapply(1:k, function(i) (1:k)[-i])))
}
make.info.mk2 <- function(phy) {
  info <- make.info.mkn(2, phy)
  info$name <- "mk2"
  info$name.pretty <- "Mk2"
  info$argnames <- default.argnames.mk2()
  info
}
default.argnames.mk2 <- function()
  c("q01", "q10")

## 3: make.cache
make.cache.mkn <- function(tree, states, k, strict, control) {
  method <- control$method

  tree <- check.tree(tree)
  if ( !is.null(states) ) # for multitrait
    states <- check.states(tree, states, strict=strict,
                           strict.vals=1:k)
  cache <- make.cache(tree)
  if ( method == "mk2" )
    cache$info <- make.info.mk2(tree)
  else
    cache$info <- make.info.mkn(k, tree)

  cache$states  <- states
  if ( method == "ode" ) {
    cache$info$derivs <- derivs.mkn.ode
    cache$y <- initial.tip.mkn.ode(cache)
    cache$info$name.ode <- "mknode"
  } else if ( method == "" ) {
    cache$info$name.ode <- "mknpij"
  }

  cache
}

## 4: initial conditions
initial.conditions.mkn <- function(init, pars, t, idx) {
  ## In a previous version I had this check here:
  ## if ( !any(ret > 0) )
  ##   stop("Incompatible initial conditions at tip ", idx)
  ## which was nice, but slows things down a bit.
  ##
  ## A standalone check might be good.
  init[,1] * init[,2]
}

## 5: rootfunc
rootfunc.mkn <- function(res, pars, root, root.p, intermediates) {
  d.root <- res$vals
  lq <- res$lq
  k <- length(d.root)

  root.p <- root.p.calc(d.root, pars, root, root.p,
                        stationary.freq.mkn)
  if ( root == ROOT.ALL )
    loglik <- log(d.root) + sum(lq)
  else
    loglik <- log(sum(root.p * d.root)) + sum(lq)

  if ( intermediates ) {
    res$root.p <- root.p
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- d.root
  }

  loglik
}

make.all.branches.mkn <- function(cache, control) {
  if ( control$method == "ode" ) {
    if ( !is.null(control$backend) && control$backend == "expokit" )
      make.all.branches.mkn.expokit(cache, control)
    else
      make.all.branches.dtlik(cache, control, initial.conditions.mkn)
  } else { # method == "pij"
    make.all.branches.mkn.pij(cache, control)
  }
}

######################################################################
## Additional functions:
stationary.freq.mkn <- function(pars) {
  ## When used in rootfunc, I think we always have a matrix here.
  stop("Needs testing")
  if ( length(pars) == 2 )
    pars[2:1] / sum(pars)
  else
    ## This is really easy (I think we just look at the eigenvalues of
    ## the Q matrix.
    .NotYetImplemented()
}

## Parameter manipulation:
## Makes a function that converts k(k-1) parameters into a k^2 Q
## matrix.
make.pars.mkn <- function(k) {
  qmat <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  npar <- k*(k-1)

  function(pars) {
    check.pars.nonnegative(pars, npar)
    qmat[idx] <- pars
    diag(qmat) <- -rowSums(qmat)
    qmat
  }
}

## This is not the most efficient possible, but should not really be
## used in code for which this is a problem.  This just wraps around
## the above function so that this can be done as a once-off.
mkn.Q <- function(pars, k) {
  if ( missing(k) )
    k <- (1 + sqrt(1 + 4*(length(pars))))/2
  if ( abs(k - round(k)) > 1e-6 || length(pars) != k*(k-1) )
    stop("Invalid parameter length")
  make.pars.mkn(k)(pars)
}

## Checking:
check.control.mkn <- function(control, k) {
  control <- modifyList(list(method="pij"), control)
  if ( control$method == "mk2" && k != 2 )
    stop("Method 'mk2' only valid when k=2")
  methods <- c("pij", "mk2", "ode")
  if ( !(control$method %in% methods) )
    stop(sprintf("control$method must be in %s",
                 paste(methods, collapse=", ")))
  control
}

check.control.mk2 <- function(control) {
  control <- modifyList(list(method="mk2"), control)
  if ( control$method != "mk2" )
    stop("Invalid control$method value (just omit it)")
  control
}
