## Currently, I have added interfaces to several different optimisers:
##   'optim' (L-BFGS-B, Nelder-Mead, BFGS, CG, SANN)
##   'subplex'
##   'nlminb'
##   'nlm'
## These now have a similar interface, though not identical.
find.mle <- function(func, x.init, method, ...) {
  UseMethod("find.mle")
}

find.mle.default <- function(func, x.init, method,
                             fail.value=NA, class.append=NULL, ...) {
  if ( is.constrained(func) )
    x.init <- guess.constrained.start(func, x.init)

  ans <- do.mle.search(func, x.init, method, fail.value, ...)
  class(ans) <- c(class.append, class(ans))
  ans
}

do.mle.search <- function(func, x.init, method, fail.value=-Inf,
                          control=list(), lower=-Inf, upper=Inf,
                          hessian=FALSE, verbose=0, keep.func=TRUE,
                          ...) {
  method <- match.arg(method, c("optim", "subplex", "nlminb", "nlm",
                                "minqa", "optimize", "int1d",
                                "mixed", "subplexR"))

  control$y.init <- y.init <- func(x.init, ...)
  if ( inherits(y.init, "try-error") )
    stop("The above error occured when testing the starting position")
  if ( !is.finite(y.init) )
    stop("Starting point must have finite probability")

  if ( is.null(control$fail.penalty) )
    control$fail.penalty <- 1000

  ## Can handle -Inf: subplex, nlminb, int1d, mixed,
  ## Require finite values: optim, nlm, minqa, optimize
  if ( is.na(fail.value) ) {
    if ( method %in% c("optim", "nlm", "minqa", "optimize") )
      fail.value <- y.init - control$fail.penalty
    else if ( method %in% c("subplex", "nlminb", "int1d", "mixed",
                            "subplexR") )
      fail.value <- -Inf
  }

  ## Protect the function, and combine in the function arguments:
  func2 <- protect(function(x) func(x, ...), fail.value)
  
  ## Add in verbosity, if needed:
  control$verbose <- verbose > 0
  if ( control$verbose )
    func2 <- func2 <- big.brother(func2, verbose)

  ## Remember the names of the input vector, or try and get it from
  ## the supplied function.
  if ( is.null(names(x.init)) && !is.null(x.init) ) {
    names.v <- try(argnames(func), silent=TRUE)
    if ( inherits(names.v, "try-error") ) {
      names.v <- NULL
    } else {
      if ( length(names.v) != length(x.init) )
        stop("Invalid parameter names length: expected ",
             length(names.v), ", got ", length(x.init))
      names(x.init) <- names.v
    }
  } else {
    names.v <- names(x.init)
  }

  mle.search <- get(sprintf("do.mle.search.%s", method))
  ans <- mle.search(func2, x.init, control, lower, upper)

  if ( verbose )
    cat("\n")

  if ( hessian ) {
    if ( !requireNamespace("numDeriv") )
      warning("The package numDeriv is required to compute the hessian")
    else
      ans$hessian <- numDeriv::hessian(func2, ans$par)
  }

  ## TODO: For constrained cases, store the full parameter vector
  ## too.
  names(ans$par) <- names.v
  ## ans$func <- func # drop or make option - this can be annoying
  ans$method <- method
  if ( is.constrained(func) ) {
    ## There used to be an option "extra" to coef.fit.mle that would
    ## do this:
    ## if ( extra && !is.null(extra.v <- attr(func, "extra")) )
    ##   c(object$par[extra.v], func(object$par, pars.only=TRUE))
    ## However, I don't know if it was ever used.
    ans$par.full <- func(ans$par, pars.only=TRUE)
  }

  ## I'm using attributes here, rather than saving as an element, for
  ## symmetry with the mcmc() function.  A better approach would
  ## probably be to define a generic function that will extract the
  ## likelihood from any fitted object.
  if (keep.func)
    attr(ans, "func") <- set.defaults(func, defaults=list(...))
  ans$func.class <- class(func) # used by anova()

  class(ans) <- "fit.mle"
  ans
}

logLik.fit.mle <- function(object, ...) {
  ll <- object$lnLik
  attr(ll, "df") <- length(object$par)
  class(ll) <- "logLik"
  ll
}

coef.fit.mle <- function(object, full=FALSE, extra=FALSE, ...) {
  if ( full && !is.null(object$par.full) )
    object$par.full
  else
    object$par
}

extractAIC.fit.mle <- function(fit, scale, k=2, ...)
  c(length(stats::coef(fit)), AIC(fit))

## Code based on MASS:::anova.negbin and ape:::anova.ace
## 
## There are two ways that I am interested in seeing models come into
## here, and within that, a couple of different possibilities:
##
## 1. A series of comparisons against a single model.  This will
##    *always* be the first model in the list.
## 
##    a. Focal model has most parameters (full) and subsequent models
##       are reduced versions nested within it.
##
##    b. Focal model has least parameters (minimal) and subsequent
##       models are generalisations of it.
##
## 2. A sequential series of comparisons (i.e., 
anova.fit.mle <- function(object, ..., sequential=FALSE) {
  mlist <- c(list(object), list(...))
  if ( length(mlist) == 1L )
    stop("Need to specify more than one model")

  cl <- sub("^fit.mle.", "", sapply(mlist, function(x) class(x)[1]))
  cl <- sub("\\..+$", "", cl)
  if ( length(unique(cl)) > 1 )
    stop("More than one class of model present -- tests invalid with mixed types")
  if ( is.null(names(mlist)) )
    names(mlist) <-
      c("", model=sprintf("model %d", seq_len(length(mlist)-1)))
  
  ## Next, extract the log-likelihoods and degrees of freedom from
  ## each model.
  ll <- lapply(mlist, stats::logLik)
  df <- sapply(ll, attr, "df")
  ll.val <- sapply(ll, as.numeric)

  ## Check if a model 'sub' appears to be nested within a model
  ## 'full'.
  check <- function(full, sub) {
    if ( all(names(stats::coef(sub)) %in% names(stats::coef(full))) ) {
      TRUE
    } else {
      drop <- c("constrained", "function")
      cl.full <- setdiff(full$func.class, drop)
      cl.sub  <- setdiff(sub$func.class,  drop)
      all(cl.sub %in% cl.full)
    }
  }

  if ( sequential ) { 
    chisq <- c(NA, 2*diff(ll.val))
    ddf <- c(NA, diff(df))
    if ( all(ddf[-1] > 0) ) {
      names(mlist)[1] <- "minimal"
      ok <- sapply(seq_len(length(mlist)-1), function(i)
                   check(mlist[[i+1]], mlist[[i]]))
    } else {
      names(mlist)[1] <- "full"
      ok <- sapply(seq_len(length(mlist)-1), function(i)
                   check(mlist[[i]], mlist[[i+1]]))
      chisq <- -chisq
      ddf <- -ddf
    }
  } else {
    chisq <- c(NA, 2*(ll.val[1] - ll.val[-1]))
    ddf <- c(NA, df[1] - df[-1])
    if ( all(df[-1] < df[1]) ) {
      ## Case a: Model 1 has the most parameters (full model) and other
      ## models should be nested within it.
      names(mlist)[1] <- "full"
      ok <- sapply(mlist[-1], check, full=mlist[[1]])
      if ( !all(ok) )
        warning(sprintf("Model(s) %s not nested within the first model?",
                        paste(which(!ok)+1, collapse=", ")))
    } else if ( all(df[-1] > df[1]) ) {
      ## Case b: Model 1 has the fewest parameters (minimal model) and
      ## is nested in all subsequent models.
      names(mlist)[1] <- "minimal"
      ok <- sapply(mlist[-1], check, sub=mlist[[1]])
      if ( !all(ok) )
        warning(sprintf("First model is not nested within model(s) %s?",
                        paste(which(!ok)+1, collapse=", ")))
      chisq <- -chisq
      ddf <- -ddf
    } else {
      stop("Incompatible models (see ?find.mle)")
    }
  }

  if ( any(chisq[-1] < 0 ) )
    warning("Impossible chi-square values: probable convergence failures")

  out <- data.frame(Df=df,
                    lnLik=sapply(ll, as.numeric),
                    AIC=sapply(mlist, AIC),
                    ChiSq=chisq,
                    "Pr(>|Chi|)"=1 - pchisq(chisq, ddf),
                    check.names=FALSE)
  rownames(out) <- names(mlist)
    
  class(out) <- c("anova", "data.frame")
  out
}

## For want of a better name, this does the initial parameter
## guessing.
guess.constrained.start <- function(func, x.init, warn=TRUE) {
  f.orig <- environment(func)$f
  names.orig <- argnames(f.orig)
  names.cons <- argnames(func)
  n.orig <- length(names.orig)
  n.cons <- length(names.cons)
  arg.idx <- match(names.cons, names.orig)

  if ( length(x.init) == n.cons ) {
    x.init
  } else if  ( length(x.init) == n.orig && !any(is.na(arg.idx)) ) {
    if ( warn )
      warning("Guessing parameters while constraining model - may do badly")
    x.init <- x.init[arg.idx]
  } else {
    stop("Could not guess reduced parameter set from those given")
  }

  x.init
}

## Simple function to drop the likelihood function from a fitted
## object.
drop.likelihood <- function(object) {
  attr(object, "func") <- NULL
  object
}

## And to retrieve it.
get.likelihood <- function(object) {
  lik <- attr(object, "func")
  if (is.null(lik) || !is.function(lik))
    stop("Did not find likelihood function in diversitree object")
  lik
}
