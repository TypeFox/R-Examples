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

## 1: make
make.musse <- function(tree, states, k, sampling.f=NULL, strict=TRUE,
                       control=list()) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions.musse)
  rootfunc <- rootfunc.musse
  f.pars <- make.pars.musse(k)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("musse", "dtlik", "function")
  ll
}

## 2: info
make.info.musse <- function(k, phy) {
  list(name="musse",        # canonical name (as in make.<name>)
       name.pretty="MuSSE", # pretty name for printing
       ## Parameters:
       np=as.integer(k * (k + 2)), # number of ODE parameters
       argnames=default.argnames.musse(k),
       ## Variables:
       ny=as.integer(2*k),            # number of variables
       k=as.integer(k),               # number of states
       idx.e=as.integer(1:k),         # index of 'E' variables
       idx.d=as.integer((k+1):(2*k)), # index of 'D' variables
       ## R version of derivative function
       derivs=derivs.musse,       
       ## Phylogeny:
       phy=phy,   # here to help with printing, possibly plotting
       ## Inference:
       ml.default="subplex", # default ML search
       mcmc.lowerzero=TRUE,  # all paramters positive? (for mcmc)
       ## These are optional
       doc=NULL,    # extra string to print during print()
       reference=c( # vector of references
         "FitzJohn (submitted)"))
}
default.argnames.musse <- function(k) {
  fmt <- sprintf("%%0%dd", ceiling(log10(k + .5)))
  str <- sprintf(fmt, 1:k)
  c(sprintf("lambda%s", str),
    sprintf("mu%s", str),
    sprintf("q%s%s", rep(str, each=k-1),
            unlist(lapply(1:k, function(i) str[-i]))))
}

## 3: make.cache (& initial.tip)
make.cache.musse <- function(tree, states, k, sampling.f=NULL,
                             strict=TRUE) {
  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=strict, strict.vals=1:k)

  cache <- make.cache(tree)
  cache$info <- make.info.musse(k, tree)
  cache$states <- states
  cache$sampling.f <- check.sampling.f(sampling.f, k)
  cache$y <- initial.tip.xxsse(cache)
  cache
}

## 4: initial.conditions:
## Note that we ignore both 't' and 'idx'.
initial.conditions.musse <- function(init, pars, t, idx) {
  k <- nrow(init)/2
  i <- seq_len(k)
  j <- i + k

  c(init[i,1],
    init[j,1] * init[j,2] * pars[i])
}

## 5: rootfunc
rootfunc.musse <- function(res, pars, condition.surv, root, root.p,
                           intermediates) {
  vals <- res$vals
  lq <- res$lq
  k <- length(vals)/2

  i <- seq_len(k)
  d.root <- vals[-i]

  ## Because this is shared with BiSSE:
  root.equi <- if ( k == 2 ) stationary.freq.bisse else NULL
  root.p <- root.p.calc(d.root, pars, root, root.p, root.equi)
  if ( condition.surv ) {
    lambda <- pars[i]
    e.root <- vals[i]
    d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2)
  }

  if ( root == ROOT.ALL )
    loglik <- log(d.root) + sum(lq)
  else
    loglik <- log(sum(root.p * d.root)) + sum(lq)

  if ( intermediates ) {
    res$root.p <- root.p
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- vals
  }

  loglik
}

###########################################################################
## Additional functions
## Heuristic starting point
starting.point.musse <- function(tree, k, q.div=5, yule=FALSE) {
  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  r <- if  ( pars.bd[1] > pars.bd[2] )
    (pars.bd[1] - pars.bd[2]) else pars.bd[1]
  p <- rep(c(pars.bd, r / q.div), c(k, k, k * (k-1)))
  names(p) <- default.argnames.musse(k)
  p
}

## For historical and debugging purposes, not used directly in the
## calculations, but branches function is generated this way
## internally.
make.branches.musse <- function(cache, control)
  make.branches.dtlik(cache$info, control)

## Parameter manipulation
make.pars.musse <- function(k) {
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.lm <- 1:(2*k)
  idx.q <- (2*k+1):(k*(1+k))

  function(pars) {
    check.pars.musse(pars, k)
    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    c(pars[idx.lm], qmat)
  }
}

## This makes the Q matrix from a set of parameters.
musse.Q <- function(pars, k) {
  if ( missing(k) )
    k <- (sqrt(1+4*length(pars))-1)/2
  if ( abs(k - round(k)) > 1e-6 || length(pars) != k*(1+k) )
    stop("Invalid parameter length")
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.q <- (2*k+1):(k*(1+k))
  qmat[idx.qmat] <- pars[idx.q]
  diag(qmat) <- -rowSums(qmat)
  qmat
}

check.pars.musse <- function(pars, k)
  check.pars.nonnegative(pars, k*(k+1))

derivs.musse <- function(t, y, pars) {
  k <- length(y)/2L
  i1 <- seq_len(k)
  i2 <- i1 + k
  i3 <- (2*k+1L):(k*(k+2L))

  E <- y[i1]
  D <- y[i2]
  lambda <- pars[i1]
  mu     <- pars[i2]
  Q      <- matrix(pars[i3], k, k)

  c(mu - (lambda + mu) * E +     lambda * E * E + Q %*% E,
       - (lambda + mu) * D + 2 * lambda * E * D + Q %*% D)
}

## Version from 9c44a294a5fd4fb55853c618d91f95eb04c79f94:
## make.musse.eqs.R <- function(k) {
##   qmat <- matrix(0, k, k)
##   idx.qmat <- cbind(rep(1:k, each=k-1),
##                unlist(lapply(1:k, function(i) (1:k)[-i])))
##   idx.e <- 1:k
##   idx.d <- (k+1):(2*k)
##   idx.l <- 1:k
##   idx.m <- (k+1):(2*k)
##   idx.q <- (2*k+1):(k*(1+k))

##   function(t, y, parms, ...) {
##     e <- y[idx.e]
##     d <- y[idx.d]
##     lambda <- parms[idx.l]
##     mu     <- parms[idx.m]
    
##     qmat[idx.qmat] <- parms[idx.q]  
##     diag(qmat) <- -rowSums(qmat)

##     list(c(mu - (lambda + mu) * e + lambda * e * e + qmat %*% e,
##            -(lambda + mu) * d + 2 * lambda * d * e + qmat %*% d))
##   }
## }
