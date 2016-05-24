initial.tip.bm.pruning <- function(cache) {
  y <- mapply(function(mean, sd) c(mean, sd*sd, 0),
              cache$states, cache$states.sd, SIMPLIFY=TRUE)
  dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
}

## 4: initial.conditions
## Note that we ignore both 't' and 'idx'.
initial.conditions.bm.pruning <- function(init, pars, t, idx) {
  m1 <- init[1,1]
  m2 <- init[1,2]
  v1 <- init[2,1]
  v2 <- init[2,2]
  vv <- v1 + v2

  c((m1 * v2 + m2 * v1) / vv,
    (v1 * v2) / vv,
    -(m1 - m2)^2 / (2 * vv) - log(2 * pi * vv) / 2)
}

## 5: rootfunc
rootfunc.bm.pruning <- function(res, pars, root, root.x,
                               intermediates) {
  vals <- res$vals
  lq <- res$lq
  
  if ( root == ROOT.MAX ) {
    ## This treats a prior on the root as a delta function centred at
    ## the ML root state.
    ## The first term can be more intuitively written as:
    ##   dnorm(vals[1], vals[1], sqrt(vals[2]), TRUE)
    ##   dnorm(0, 0, sqrt(vals[2]), TRUE)
    loglik <- - log(2 * pi * vals[[2]]) / 2 + vals[[3]] + sum(lq)
  } else if ( root == ROOT.FLAT ) {
    ## Flat prior (by this point, function integrates to vals[[3]])
    loglik <- vals[[3]] + sum(lq)
  } else if ( root == ROOT.OBS ) {
    ## Observed weighting (integrate norm norm wrt x from -inf to inf
    ## gives 1 / (2 sqrt(pi s2))).
    loglik <- -log(2 * sqrt(pi * vals[[2]])) + vals[[3]] + sum(lq)
  } else if ( root == ROOT.GIVEN ) {
    if ( is.null(root.x) )
      stop("root.x not provided, but root=ROOT.GIVEN specified")
    loglik <- dnorm(root.x, vals[1], sqrt(vals[2]), TRUE) +
      vals[[3]] + sum(lq)
  } else {
    stop("Invalid root mode")
  }

  if ( intermediates ) {
    res$root.p <- NA # not sure what would be good here...
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- vals
  }

  loglik
}

make.all.branches.bm.pruning <- function(cache, control) {
  if ( control$backend == "R" )
    function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache,
                          initial.conditions.bm.pruning,
                          branches.bm.pruning, preset)
  else # backend == "C"
    make.all.branches.continuous(cache, control)
}

###########################################################################
## Additional functions

## branches
## Unlike the ODE-based functions, this ignores t0 and carries out
## calculations along length=len (rather than looking at c(t0,
## t0+len)).
branches.bm.pruning <- function(y, len, pars, t0, idx) {
  m <- y[1]
  v <- y[2]
  z <- y[3]
  list(z, c(m, v + (pars[1] * len), 0))
}

