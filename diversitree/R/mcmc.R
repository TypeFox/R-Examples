## Slice sampling

## These functions take a function 'f' that evaluates an argument 'x'
## which is a vector of parameters.  Across a single iteration,
## take the input x and return a vector of parameters 'y'
## corresponding to a new position.

## The default MCMC method will return a half-finished MCMC sample if
## interrupted with Ctrl-C.  It is expected that the default method
## will be sufficiently general to be useful for most approaches.
## Non-default methods should normally be default parameter setting.
## However, totally different MCMC algorithms could be used.

## Currently, the 'control' argument is not used by any method.
mcmc <- function(lik, x.init, nsteps, ...) {
  UseMethod("mcmc")
}

## This function is starting to get quite unwieldly; I'll probably
## split it into a preparation and main section at some point.
mcmc.default <- function(lik, x.init, nsteps, w, prior=NULL,
                         sampler=sampler.slice, fail.value=-Inf,
                         lower=-Inf, upper=Inf, print.every=1,
                         control=list(),
                         save.file, save.every=0, save.every.dt=NULL,
                         previous=NULL, previous.tol=1e-4,
                         keep.func=TRUE,
                         ...) {
  if ( is.null(sampler) )
    sampler <- sampler.slice
  
  if ( save.every > 0 || !is.null(save.every.dt) ) {
    if ( missing(save.file) )
      stop("save.file must be given if save.every > 0")
    save.type <- tools::file_ext(save.file)
    save.type <- match.arg(tolower(save.type), c("csv", "rds"))
    save.fun <- switch(save.type,
                       rds=saveRDS,
                       csv=function(x, f) write.csv(x, f, row.names=FALSE))
    save.file.bak <- paste(save.file, ".bak", sep="")
  }

  n.previous <- if ( is.null(previous) ) 0 else nrow(previous)
  if ( !is.null(previous) ) {
    if ( !inherits(previous, "mcmcsamples") )
      stop("Currently only mcmcsamples objects can be continued")
    if ( n.previous >= nsteps ) {
      warning("Chain already complete")
      return(previous)
    }
    if ( !is.null(x.init) )
      stop("x.init must be NULL if continuing")
    npar <- ncol(stats::coef(previous)) # tons of repetition here.
    hist.pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist.prob <- rep(NA, nsteps)
    
    hist.pars[seq_len(n.previous),] <- stats::coef(previous)
    hist.prob[seq_len(n.previous)]  <- previous$p

    x.init <- hist.pars[n.previous,]
    y.prev <- hist.prob[n.previous]
  } else {
    npar <- length(x.init)
    hist.pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist.prob <- rep(NA, nsteps)
  }

  if ( is.null(names(x.init)) )
    try(colnames(hist.pars) <- names(x.init) <- argnames(lik),
        silent=TRUE)
  else
    colnames(hist.pars) <- names(x.init)
  
  if ( is.null(prior) )
    posterior <- protect(function(x) lik(x, ...),
                         fail.value.default=fail.value)
  else
    posterior <- protect(function(x) lik(x, ...) + prior(x),
                         fail.value.default=fail.value)

  lower <- check.par.length(lower, npar)
  upper <- check.par.length(upper, npar)
  w     <- check.par.length(w,     npar)

  check.bounds(lower, upper, x.init)

  y.init <- posterior(x.init, fail.value=NULL)
  if ( !is.finite(y.init) ||
      (!is.null(fail.value) && y.init == fail.value) )
    stop("Starting point must have finite probability")
  if ( n.previous > 0 && abs(y.prev - y.init) > previous.tol ) {
    msg <- paste("Cannot continue the chain:\n",
                 sprintf("\texpected posterior = %2.7, got %2.7f",
                         previous$p[n.previous], y.init))
    stop(msg)
  }

  class.str <- c(sprintf("mcmcsamples.%s", get.info(lik)$name),
                 "mcmcsamples", "data.frame")

  clean.hist <- function(pars, p) {
    out <- data.frame(i=seq_along(p), pars, p)
    class(out) <- class.str
    out
  }

  we.should.print <- make.every.so.often(print.every)
  we.should.save <- make.every.so.often(save.every, save.every.dt)

  mcmc.loop <- function() {
    for ( i in seq(n.previous+1, nsteps, by=1) ) {
      tmp <- sampler(posterior, x.init, y.init, w,
                     lower, upper, control)
      x.init <- hist.pars[i,] <<- tmp[[1]]
      y.init <- hist.prob[i]  <<- tmp[[2]]

      if ( we.should.print() )
        cat(sprintf("%d: {%s} -> %2.5f\n", i,
                    paste(sprintf("%2.4f", x.init), collapse=", "),
                    y.init))
      if ( we.should.save() ) {
        j <- seq_len(i)
        ## Back up the old version to avoid IO errors if the system
        ## fails while saving.
        if ( file.exists(save.file) )
          ok <- file.rename(save.file, save.file.bak)
        ok <- try(save.fun(clean.hist(hist.pars[j,,drop=FALSE], hist.prob[j]),
                           save.file))
        if ( inherits(ok, "try-error") )
          warning("Error while writing progress file (continuing)",
                  immediate.=TRUE)
      }
    }
    clean.hist(hist.pars, hist.prob)
  }

  mcmc.recover <- function(...) {
    j <- !is.na(hist.prob)
    if ( !any(j) )
      stop("MCMC was stopped before any samples taken")
    hist <- clean.hist(hist.pars[j,], hist.prob[j])
    warning("MCMC was stopped prematurely: ", nrow(hist), "/", nsteps,
            " steps completed.  Truncated chain is being returned.",
            immediate.=TRUE)
    hist
  }

  samples <- tryCatch(mcmc.loop(), interrupt=mcmc.recover)

  if ( save.every > 0 || !is.null(save.every.dt) )
    if ( nrow(samples) == nsteps && file.exists(save.file.bak) )
      file.remove(save.file.bak)

  if (keep.func) {
    attr(samples, "func")  <- set.defaults(lik, defaults=list(...))
    attr(samples, "prior") <- prior
  }

  samples
}

mcmc.dtlik <- function(lik, x.init, nsteps, lower=-Inf, ...) {
  if ( missing(lower) && get.info(lik)$mcmc.lowerzero )
    lower <- 0
  NextMethod("mcmc", lower=lower)
}

## This is common, so this helps reduce code duplication.
mcmc.lowerzero <- function(lik, x.init, nsteps, ..., lower=0)
  NextMethod("mcmc", lower=lower)

make.unipar <- function(f, x, i) {
  function(z) {
    x[i] <- z
    f(x)
  }
}

make.prior.exponential <- function(r) {
  function(pars)
    sum(dexp(pars, r, log=TRUE))
}

## This is still experimental, and will not work nicely unless
## everything is nicely paired (it will not work well with constrained
## models, for example).
make.prior.ExpBeta <- function(r, beta) {
  to.pars2 <- function(pars) {
    m <- matrix(pars, 2)
    pars.mean <- colMeans(m)
    d <- 1 - (m[1,] / (2*pars.mean))
    rbind(pars.mean, d)
  }
  function(pars) {
    pars2 <- to.pars2(pars)
    sum(dexp(pars2[1,], r, log=TRUE)) +
      sum(dbeta(pars2[2,], beta, beta, log=TRUE))
  }
}

## TODO: Allow vector lower and upper here...
make.prior.uniform <- function(lower, upper, log=TRUE) {
  if ( length(lower) == 2 && missing(upper) ) {
    upper <- lower[2]
    lower <- lower[1]
  }
  n <- length(lower)
  if ( length(upper) != n )
    stop("'lower' and 'upper' both be the same length")
  p.in <- 1/(upper - lower)
  p.out <- 0
  if ( log ) {
    p.in <- log(p.in)
    p.out <- -Inf
  }
  function(x) {
    ret <- rep(p.in, length.out=length(x))
    ret[x < lower | x > upper] <- p.out
    if ( log )
      sum(ret)
    else
      prod(ret)
  }
}

mcmcsamples.index <- function(n, burnin=NA, thin=NA, sample=NA) {
  i <- seq_len(n)
  if (!is.na(burnin) && burnin > 0) {
    if (burnin < 1)
      burnin <- floor(burnin * n)
    i <- i[i > burnin]
  }
  if (!is.na(thin) && thin > 1)
    i <- i[seq(1, length(i), by=thin)]
  if (!is.na(sample)) {
    if (sample > length(i))
      warning("Sampling *will* generate duplicates")
    i <- i[sample(length(i), sample, replace=TRUE)]
  }
  i
}

coef.mcmcsamples <- function(object, burnin=NA, thin=NA, sample=NA,
                             full=FALSE, ...) {
  i <- mcmcsamples.index(nrow(object), burnin, thin, sample)
  p <- as.matrix(object[i,-c(1, ncol(object)),drop=FALSE])
  if (full) {
    lik <- get.likelihood(object)
    if (inherits(lik, "constrained"))
      p <- t(apply(p, 1, lik, pars.only=TRUE))
  }
  p
}

as.mcmcsamples <- function(x, ...) {
  if ( !identical(colnames(x)[c(1, ncol(x))], c("i", "p")) )
    stop("Nope")
  x <- as.data.frame(x)
    class(x) <- c("mcmcsamples", "data.frame")
  x
}

## This is a bit of an overkill because
make.every.so.often <- function(iterations=1, dt=NULL) {
  i <- 0 # counter, will be updated
  
  if ( !is.null(dt) ) {
    requireNamespace("lubridate")
    if ( !inherits(dt, "Period") )
      stop("dt must be a Period object")
    t.next <- lubridate::now() + dt
    ## Note use of 'or' here, not 'and'.
    check <- function() {
      i <<- i + 1
      ok <- (iterations > 0 && i %% iterations == 0) || lubridate::now() > t.next
      if ( ok )
        t.next <<- lubridate::now() + dt
      ok
    }
  } else if ( iterations > 0 ) {
    check <- function() {
      i <<- i + 1
      i %% iterations == 0
    }
  } else {
    check <- function()
      FALSE
  }
  check
}

