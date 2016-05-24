## Nee 1994 equations.
make.cache.bd.nee <- function(tree=NULL, sampling.f=NULL,
                              unresolved=NULL, times=NULL) {
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  sampling.f <- check.sampling.f(sampling.f, 1)
  unresolved <- check.unresolved.bd(tree, unresolved)

  if ( is.null(times) ) { # proper tree
    tree <- check.tree(tree)
    times <- branching.times(tree)
  } else { ## Special case...
    if ( !is.null(tree) )
      stop("Can only supply times if tree=NULL")
  }

  times <- as.numeric(sort(times, decreasing=TRUE))
  x <- c(NA, times)
  N <- length(x)
  n.node <- length(times) # == tree$Nnode
  const <- lfactorial(N - 1)

  list(N=N, x=x, f=sampling.f, unresolved=unresolved, n.node=n.node,
       const=const)
}

make.all.branches.bd.nee <- function(cache, control) {
  N <- cache$N
  x <- cache$x
  f <- cache$f
  unresolved <- cache$unresolved

  sum.x.3.N <- sum(x[3:N])

  function(pars, intermediates, preset=NULL) { # preset ignored 
    if ( pars[2] == pars[1] )
      pars[1] <- pars[1] + 1e-12
    r <- pars[[1]] - pars[[2]]
    a <- pars[[2]] / pars[[1]]

    if ( f < 1 )
      loglik <-
        (N - 2) * log(f*abs(r)) + N * log(abs(1 - a)) + r * sum.x.3.N -
          2*sum(log(abs(f*exp(r * x[2:N]) - a + 1 - f)))
    else
      loglik <-
        (N - 2) * log(  abs(r)) + N * log(abs(1 - a)) + r * sum.x.3.N -
          2*sum(log(abs(  exp(r * x[2:N]) - a        )))
    if ( !is.null(unresolved) ) {
      ert <- exp(r * unresolved$t)
      loglik <- loglik +
        sum((unresolved$n-1) * (log(abs(ert - 1)) - log(abs(ert - a))))
    }

    ## At some point, this should be e.root only.  See derivations.pdf
    e <- log(f * f * r * (1 - a)) -
      2*log(abs(exp(-r * x[2])*(a - 1 + f) - f))
    list(vals=c(loglik, e))
  }
}

rootfunc.bd.nee <- function(res, pars, condition.surv,
                            intermediates, const) {
  if ( intermediates )
    stop('intermediates cannot be produced -- use method="ode"')
  vals <- res$vals
  loglik <- vals[[1]] + const
  if ( !condition.surv ) # notice this is opposite to usual!
    loglik <- loglik + vals[[2]]
  loglik
}
