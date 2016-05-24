## ODE interface to BD models.
make.cache.bd.ode <- function(tree, sampling.f, unresolved) {
  if ( is.null(tree) )
    stop('Can only supply times if method="nee"')

  tree <- check.tree(tree)
  unresolved <- check.unresolved.bd(tree, unresolved)  

  ## This is because it's not obvious what an ODE solution to the
  ## unresolved clades would be...
  if ( !is.null(unresolved) )
    stop("Cannot deal with unresolved clades yet with this method.")
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, 1)

  cache <- make.cache(tree)
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f
  cache$y <- initial.tip.bd.ode(cache)
  cache$info <- make.info.bd(tree)
  cache$const <- lfactorial((length(cache$len) + 1)/2 - 1)
  cache
}
initial.tip.bd.ode <- function(cache) {
  f <- cache$sampling.f
  y <- list(c(1-f, f))
  y.i <- rep(1, length(cache$tips))
  dt.tips.grouped(y, y.i, cache)
}

## 4: initial.conditions
## Note that we ignore both 't' and 'idx'.
initial.conditions.bd.ode <- function(init, pars, t, idx)
  c(init[1,1],
    init[2,1] * init[2,2] * pars[1])

rootfunc.bd.ode <- function(res, pars, condition.surv,
                            intermediates, const) {
  vals <- res$vals
  lq <- res$lq
  d.root <- vals[2]

  ## Compute N! for comparability with the non-ode method
  ## const <- lfactorial((length(lq) + 1)/2 - 1)
  
  if ( condition.surv ) {
    e.root <- vals[[1]]
    lambda <- pars[[1]]
    d.root <- d.root/(lambda * (1 - e.root)^2)
  }
  loglik <- log(d.root) + sum(lq) + const
  names(loglik) <- NULL

  if ( intermediates ) {
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- vals
  }

  loglik
}

######################################################################
## Additional functions:
make.branches.bd.ode <- function(cache, control)
  make.branches.dtlik(cache$info, control)  

derivs.bd <- function(t, y, pars) {
  E <- y[1]
  D <- y[2]
  lambda <- pars[1]
  mu <- pars[2]
  c(mu - (mu + lambda)*E +   lambda*E*E,
    - (mu + lambda)*D + 2*lambda*D*E)
}
