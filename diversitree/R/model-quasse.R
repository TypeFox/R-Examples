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
make.quasse <- function(tree, states, states.sd, lambda, mu,
                            control=NULL, sampling.f=NULL) {
  cache <- make.cache.quasse(tree, states, states.sd, lambda, mu,
                             control, sampling.f)
  all.branches <- make.all.branches.quasse(cache, cache$control)
  rootfunc <- make.rootfunc.quasse(cache)
  f.pars <- make.pars.quasse(cache)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.f=NULL, intermediates=FALSE) {
    pars2 <- f.pars(pars)
    ans <- all.branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.f, intermediates)
  }

  class(ll) <- c("quasse", "dtlik", "function")
  ll
}

## 2: info
make.info.quasse <- function(lambda, mu, phy) {
  ## Work around for .split:
  if ( !is.null(lambda) )
    argnames <- default.argnames.quasse(lambda, mu)
  else
    argnames <- NULL
  list(name="quasse",
       name.pretty="QuaSSE",
       ## Parameters:
       np=NA,
       argnames=argnames,
       ## Variables:
       ny=NA,
       k=NA,
       idx.e=NA,
       idx.d=NA,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=FALSE, # not for many models.
       ## These are optional
       doc=NULL,
       reference=c("FitzJohn (2010) doi:10.1093/sysbio/syq053"),
       ## These are special to QuaSSE:
       lambda=lambda,
       mu=mu)
}
default.argnames.quasse <- function(lambda, mu) {
  c(sprintf("l.%s", names(formals(lambda))[-1]),
    sprintf("m.%s", names(formals(mu))[-1]),
    "drift", "diffusion")
}

## 3: make.cache (& initial conditions)
make.cache.quasse <- function(tree, states, states.sd, lambda, mu,
                              control, sampling.f, for.split=FALSE) {
  ## 1: tree
  tree <- check.tree(tree)  

  ## 2: states & errors
  tmp <- check.states.quasse(tree, states, states.sd)
  states <- tmp$states
  states.sd <- tmp$states.sd

  ## 3: Control structure (lots of checking!)
  control <- check.control.quasse(control, tree, states)

  cache <- make.cache(tree)
  cache$states  <- states
  cache$states.sd <- states.sd
  cache$control <- control

  ## This is a bit ugly, but only do these checks if we are *not*
  ## doing a split QuaSSE model.  Function checking is done separately
  ## there, but everything above is the same.
  if ( !for.split ) {
    ## 4: Speciation/extinction functions
    n.lambda <- check.f.quasse(lambda)
    n.mu     <- check.f.quasse(mu)
    n.args   <- n.lambda + n.mu + 2
    args <- list(lambda=seq_len(n.lambda),
                 mu=seq_len(n.mu) + n.lambda,
                 drift=n.lambda + n.mu + 1,
                 diffusion=n.lambda + n.mu + 2)

    cache$lambda <- lambda
    cache$mu <- mu
    cache$args <- args

    sampling.f <- check.sampling.f(sampling.f, 1)
    cache$sampling.f <- sampling.f
  }
  cache$info <- make.info.quasse(lambda, mu, tree)

  cache
}
initial.tip.quasse <- function(cache, control, x) {
  nx <- control$nx * control$r
  npad <- nx - length(x)
  e0 <- 1 - cache$sampling.f

  if ( control$tips.combined ) {
    tips <- cache$tips
    t <- cache$len[tips]
    i <- order(t)
    target <- tips[i]

    states <- cache$states[i]
    states.sd <- cache$states.sd[i]

    y <- mapply(function(mean, sd)
                c(dnorm(x, mean, sd), rep(0, npad)),
                states, states.sd, SIMPLIFY=FALSE)
    y <- matrix(c(rep(e0, nx), unlist(y)), nx, length(target)+1)

    list(target=target, y=y, t=t[i])
  } else {
    y <- mapply(function(mean, sd)
                c(rep(e0, nx),
                  dnorm(x, mean, sd), rep(0, npad)),
                cache$states, cache$states.sd, SIMPLIFY=FALSE)
    dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
  }
}

## 4: initial.conditions
make.initial.conditions.quasse <- function(control) {
  tc <- control$tc
  r <- control$r
  nx.lo <- control$nx
  nx.hi <- nx.lo * r

  ## There is the chance that we could be slightly off on the depth
  ## by rounding error.  Because of this, I've done the testing
  ## against the *length* of the data, and then checked that the time
  ## is appropriate (to within eps of the correct value).  It is
  ## possible that two different branches with different numbers of
  ## nodes that finish right at the critical interval might have
  ## conflicting lengths.
  eps <- 1e-8
  function(init, pars, t, idx) {
    if ( length(init[[1]]) != length(init[[2]]) )
      stop("Data have incompatible length")

    if ( t < tc ) {
      ## if ( length(init[[1]]) / 2 == nx.hi ) { # t < tc
      ## if ( !((t - eps) < tc) )
      ##   stop("Wrong data size")
      nx <- nx.hi
      lambda <- pars[[1]]$lambda
    } else {
      ## if ( !((t + eps) > tc) )
      ##   stop("Wrong data size")
      nx <- nx.lo
      lambda <- pars[[2]]$lambda
    }
    
    ndat <- length(lambda)
    i <- seq_len(nx)
    j <- seq.int(nx+1, nx + ndat)

    c(init[[1]][i],
      init[[1]][j] * init[[2]][j] * lambda,
      rep.int(0.0, nx - ndat))
  }
}

## 5 rootfunc
## This function assumes that the root node is in the low-condition,
## which is enforced by the checking.
make.rootfunc.quasse <- function(cache) {
  root.idx <- cache$root
  nx <- cache$control$nx
  dx <- cache$control$dx
  
  function(res, pars, condition.surv, root, root.f, intermediates) {
    vals <- matrix(res$vals, nx, 2)[seq_len(pars$lo$ndat),]
    lq <- res$lq

    d.root <- vals[,2]

    root.p <- root.p.quasse(d.root, pars$lo, root, root.f)
    if ( condition.surv ) {
      lambda <- pars$lo$lambda
      e.root <- vals[,1]
      d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2) * dx
    }

    loglik <- log(sum(root.p * d.root) * dx) + sum(lq)

    if ( intermediates ) {
      attr(loglik, "intermediates") <- res
      attr(loglik, "vals") <- vals
    }
    loglik
  }
}

root.p.quasse <- function(d.root, pars, root, root.f) {
  if ( !is.null(root.f) && root != ROOT.GIVEN )
    warning("Ignoring specified root state")
  
  x <- pars$x
  dx <- x[2] - x[1]

  if ( root == ROOT.FLAT ) {
    p <- 1 / ((pars$nx-1) * dx)
  } else if ( root == ROOT.OBS )  {
    p <- d.root / (sum(d.root) * dx)
  } else if ( root == ROOT.GIVEN ) {
    p <- root.f(x)
  } else {
    stop("Unsupported root mode")
  }

  p
}

######################################################################
## Extra core stuff:
make.all.branches.quasse <- function(cache, control) {
  branches <- make.branches.quasse(cache, control)
  initial.conditions <- make.initial.conditions.quasse(control)
  ## TODO: This is where tips.combined goes, *not* in the likelihood
  ## function...
  function(pars, intermediates, preset=NULL) {
    cache$y <- initial.tip.quasse(cache, cache$control, pars[[1]]$x)
    all.branches.list(pars, cache, initial.conditions,
                      branches, preset)
  }
}

make.branches.quasse <- function(cache, control) {
  ## TODO: strictly, these should be backends...
  if ( control$method == "fftC" )
    branches <- make.branches.quasse.fftC(control)
  else if ( control$method == "fftR" )
    branches <- make.branches.quasse.fftR(control)
  else if ( control$method == "mol" )
    branches <- make.branches.quasse.mol(control)
  else # already checked.
    stop("Unknown method", control$method)
}

check.pars.quasse <- function(lambda.x, mu.x, drift, diffusion) {
  if ( any(!is.finite(c(lambda.x, mu.x, drift, diffusion))) )
    stop("Non-finite/NA parameters")
  if ( any(lambda.x < 0) || any(mu.x < 0) || diffusion <= 0 )
    stop("Illegal negative parameters")
  if ( !any(lambda.x > 0) )
    stop("No positive lambda; cannot compute likelihood")
}


## Huge chunks of this are shared with predict.dtlik.t, but it's not
## clear yet where the similarities lie.
## predict.quasse <- function(object, p, x, nx=101, v=NULL,
##                            thin=10, alpha=1/20, ...) {
##   cache <- get.cache(object)

##   if ( inherits(p, "fit.mle") )
##     p <- stats::coef(p, full=TRUE)
##   else if ( inherits(p, "mcmcsamples") )
##     p <- stats::coef(p, lik=object, full=TRUE, thin=thin)
##   ## The other case to deal with here would be constrained functions
##   ## where parameters still need expanding...
  
##   if ( missing(x) ) {
##     r <- range(cache$states, na.rm=TRUE)
##     x <- seq(r[1], r[2], length.out=nx)
##   }
##   if ( is.null(v) )
##     v <- c("lambda", "mu")

##   f <- function(i) {
##     g <- function(x, p, ...)
##       do.call(cache$info[[i]], as.list(x=x, p))
##     average.over.mcmc(p[,cache$args[[i]]], g,
##                       x, alpha=alpha)
##   }

##   for ( i in v ) {
    
##   }

##   stop()
##   ## average.over.mcmc(p, 

##   f <- function(x) {
##     y <- sort(linear(x, samples$l.c, samples$l.m))
##     c(mean(y),                            # mean
##       quantile(y, c(alpha/2, 1-alpha/2))) # range
##   }
##   list(x=xx, y=do.call(rbind, lapply(xx, f)))
## }

## predict.quasse.split <- function(object, p, x, nx, ...) {
##   stop()
## }
