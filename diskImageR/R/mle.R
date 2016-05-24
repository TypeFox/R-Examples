## These functions are from the package diversitree
## As of April 3, 2015, I have been unable to load diversitree from CRAN or github
## So, attemping to load just this set of functions, which are required.
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
  # if ( control$verbose )
    # func2 <- func2 <- big.brother(func2, verbose)

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

  # if ( hessian ) {
    # if ( !require(numDeriv) )
      # warning("The package numDeriv is required to compute the hessian")
    # else
      # ans$hessian <- hessian(func2, ans$par)
  # }

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
  c(length(coef(fit)), AIC(fit))

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
  ll <- lapply(mlist, logLik)
  df <- sapply(ll, attr, "df")
  ll.val <- sapply(ll, as.numeric)

  ## Check if a model 'sub' appears to be nested within a model
  ## 'full'.
  check <- function(full, sub) {
    if ( all(names(coef(sub)) %in% names(coef(full))) ) {
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

argnames <- function(x, ...)
  UseMethod("argnames")
`argnames<-` <- function(x, value)
  UseMethod("argnames<-")

argnames.constrained <- function(x, ...)
  attr(x, "argnames")
`argnames<-.constrained` <- function(x, value)
  stop("Cannot set argnames on constrained function")

is.constrained <- function(x)
  inherits(x, "constrained")

## First up, consider the one-shot case: don't worry about incremental
## updates.
## 
## For the first case, everything is OK on the lhs and rhs
## For subsequent cases:
##   lhs cannot contain things that are
##      - constrained things (already lhs anywhere)
##      - things constrained to (things on the rhs anywhere)
##   rhs cannot contain things that are
##      - constrained things (already lhs anywhere)
## It is possibly worth pulling out all the numerical constants and
## the "paired" parameters here to avoid using eval where
## unnecessary.  However, this makes the function substantially uglier
## for a very minor speedup.
constrain <- function(f, ..., formulae=NULL, names=argnames(f),
                      extra=NULL) {
  if ( is.constrained(f) ) {
    formulae <- c(attr(f, "formulae"), formulae)
    f <- attr(f, "func")
  }

  formulae <- c(formulae, list(...))
  names.lhs <- names.rhs <- names
  rels <- list()
  
  for ( formula in formulae ) {
    res <- constrain.parse(formula, names.lhs, names.rhs, extra)
    if ( attr(res, "lhs.is.target") ) {
      i <- try(which(sapply(rels, function(x) identical(x, res[[1]]))),
               silent=TRUE)
      if ( inherits(i, "try-error") )
        stop(sprintf("Error parsing constraint with %s on lhs",
                     as.character(res[[1]])))
      rels[i] <- res[[2]]

      ## This will not work with *expressions* involving the LHS; that
      ## would require rewriting the expressions themselves (which
      ## would not be too hard to do).  But for now let's just cause
      ## an error...
      lhs.txt <- as.character(res[[1]])
      if ( any(sapply(rels, function(x) lhs.txt %in% all.vars(x))) )
        stop(sprintf("lhs (%s) is in an expression and can't be constrained",
                     lhs.txt))
    }
    
    names.lhs <- setdiff(names.lhs, unlist(lapply(res, all.vars)))
    names.rhs <- setdiff(names.rhs, as.character(res[[1]]))
    rels <- c(rels, structure(res[2], names=as.character(res[[1]])))
  }

  i <- match(unique(sapply(rels, as.character)), extra)
  final <- c(extra[sort(i[!is.na(i)])], names.rhs)
  npar <- length(final)

  ## "free" are the parameters that have nothing special on their RHS
  ## and are therefore passed directly through
  free <- setdiff(names.rhs, names(rels))
  free.i <- match(free, names) # index in full variables
  free.j <- match(free, final) # index in given variables.

  ## Targets are processed in the same order as given by formulae. 
  target.i <- match(names(rels), names)

  pars.out <- rep(NA, length(names))
  names(pars.out) <- names
  g <- function(pars, ..., pars.only=FALSE) {
    if ( length(pars) != npar )
      stop(sprintf("Incorrect parameter length: expected %d, got %d",
                   npar, length(pars)))

    pars.out[free.i] <- pars[free.j]
    e <- structure(as.list(pars), names=final)
    pars.out[target.i] <- unlist(lapply(rels, eval, e))

    if ( pars.only )
      pars.out
    else
      f(pars.out, ...)
  }

  class(g) <- c("constrained", class(f))
  attr(g, "argnames") <- final
  attr(g, "formulae") <- formulae
  attr(g, "extra") <- extra
  attr(g, "func") <- f
  g
}


## Parsing constraints:
## The LHS of a formula must be a single variable name that exists in
## "names.lhs"
##
## The RHS can be one of
##   - numeric value
##   - expression
## If it is an expression, then all variable names must be found in
## names.rhs (or perhaps in the containing environment - check in the
## future?)

## I might eventually allow formulae of the form
##   lambda | lambda1 ~ lambda0
## to allow renaming?

## How do I use this to allow recasting to alternative bases?
## Forward definitions for just the diversification rate are
##   foo(f, r0 ~ lambda0 - mu0, r1 ~ lambda1 - mu1)
## and for both diversification and relative extinction    
##   foo(f,
##       r0 ~ lambda0 - mu0, r1 ~ lambda1 - mu1,
##       e0 ~ mu0 / lambda0, e1 ~ mu1 / lambda1)
## Backward definitions (leaving mu unchanged)
##   foo(f, lambda0 ~ r0 + mu0, lambda1 ~ r1 + mu1)
## both:
##   foo(f, lambda0 ~ r0/(1 - e0), lambda1 ~ r1/(1 - e1),
##       r0 * e0 / (1 - e0), r1 * e1 / (1 - e1))
constrain.parse <- function(formula, names.lhs, names.rhs,
                            extra=NULL) {
  formula <- as.formula(formula)
  if ( length(formula) != 3L )
    stop("Invalid formula")
  lhs <- formula[[2]]
  rhs <- formula[[3]]

  ## Checking the lhs is easy: is the lhs in the list of allowable
  ## names and of length 1?  Anything that does not match this is
  ## invalid.
  if ( !is.name(lhs) )
    stop("Invalid target on LHS of formula" )
  lhs.is.target <- is.na(match(as.character(lhs), names.lhs))

  ## Checking the rhs is more difficult.  We are OK if any of the
  ## following is met:
  ##   Numeric values (checked at the end)
  ##   If all vars in rhs are in names.rhs
  ##   There is a single var and it is in extra
  ## Failing that, if the rhs is a single variable that does exist in
  ## the calling environment.
  if ( is.language(rhs) ) {
    vars <- all.vars(rhs)
    ok <- (all(vars %in% names.rhs) ||
           length(vars) == 1 && vars %in% extra)
    if ( !ok && length(vars) == 1 ) {
      e <- parent.frame()
      if ( exists(vars, e) ) {
        rhs <- get(vars, e)
        ok <- TRUE
      }
    }

    if ( !ok )
      stop("Invalid RHS of formula:\n\t", as.character(rhs))
    if ( as.character(lhs) %in% vars )
      stop("LHS cannot appear in RHS")
  } else if ( !is.numeric(rhs) ) {
    stop("RHS must be expression, variable or number")
  }
  res <- list(lhs, rhs)
  attr(res, "lhs.is.target") <- lhs.is.target
  res
}

constrain.i <- function(f, p, i.free) {
  npar <- length(i.free)
  argnames <- try(argnames(f), silent=TRUE)
  if ( inherits(argnames, "try-error") )
    argnames <- NULL
  else if ( length(p) != length(argnames) )
    stop(sprintf("Incorrect parameter length: expected %d, got %d",
                 length(argnames), length(p)))

  pars.out <- p
  g <- function(pars, ..., pars.only=FALSE) {
    if ( length(pars) != npar )
      stop(sprintf("Incorrect parameter length: expected %d, got %d",
                   npar, length(pars)))
    pars.out[i.free] <- pars
    if ( pars.only )
      pars.out
    else
      f(pars.out, ...)
  }

  class(g) <- c("constrained.i", "constrained", class(f))
  if ( !is.null(argnames) )
    attr(g, "argnames") <- argnames[i.free]
  attr(g, "func") <- f
  g
}

## Take a function 'trans' that converts from one parameter vector to
## another and make a constrained version of a likelihood function
## 'lik' that basically evalulates
##   lik(trans(pars))
## 'argnames' must contain the names of the parameters of the
## constrained function.
do.constrain <- function(lik, trans, argnames) {
  if ( inherits(lik, "constrained") )
    stop("Cannot use do.constrain() with a constrained function")
  
  ret <- function(pars, ..., pars.only=FALSE) {
    if ( pars.only )
      trans(pars)
    else
      lik(trans(pars), ...)
  }

  class(ret) <- c("constrained", class(lik))
  attr(ret, "argnames") <- argnames
  attr(ret, "func") <- lik
  ret
}

constrain.i2 <- function(f, p, i.free) {
  npar <- length(i.free)
  argnames <- try(argnames(f), silent=TRUE)
  if ( inherits(argnames, "try-error") ) {
    argnames.constrained <- argnames <- NULL
  } else {
    argnames.constrained <- argnames[i.free]
  
    if ( length(p) != length(argnames) )
      stop(sprintf("Incorrect parameter length: expected %d, got %d",
                   length(argnames), length(p)))
    pars.out <- p
  }

  trans <- function(pars) {
    if ( length(pars) != npar )
      stop(sprintf("Incorrect parameter length: expected %d, got %d",
                   npar, length(pars)))
    pars.out[i.free] <- pars
    pars.out
  }

  ret <- do.constrain(f, trans, argnames[i.free])
  ## class(ret) <- c("constrained.i", class(ret))
  ret
}

## This is a nice idea, but some AI around base parameters would be
## nice.  In particular if there is some pattern '(X).(Y)' where X is
## unique within p?
expand.parameters <- function(p, lik.new, repl=0,
                              target=argnames(lik.new)) {
  if ( !all(names(p) %in% target) )
    stop("Not all parameters in 'p' are present in lik.new()")
  p.new <- p[target]
  names(p.new) <- target
  i <- setdiff(target, names(p.new))
  if ( length(repl) > 1 )
    if ( length(i) != length(repl) )
      stop("Wrong number of replacement parameters")
  p.new[setdiff(target, names(p))] <- repl
  p.new
}

protect <- function(f, fail.value.default=NULL) {
  function(..., fail.value=fail.value.default, finite=TRUE) {
    if ( is.null(fail.value) )
      f(...)
    else {
      ret <- try(f(...), silent=TRUE)
      failed <- (inherits(ret, "try-error") ||
                 (finite && !is.finite(ret)))
      if ( failed )
        fail.value
      else
        ret
    }
  }
}

do.mle.search.optim <- function(func, x.init, control, lower, upper) {
  control <- modifyList(list(fnscale=-1,
                             ndeps=rep(1e-5, length(x.init)),
                             optim.method="L-BFGS-B"), control)

  optim.method <- control$optim.method
  control.optim <- control[c("fnscale", "ndeps")]

  ans <- optim(x.init, func, method=optim.method,
               control=control.optim, lower=lower, upper=upper)
  names(ans)[names(ans) == "value"] <- "lnLik"  
  ans$optim.method <- optim.method
  
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (optim): ",
            tolower(ans$message))

  ans
}

do.mle.search.subplex <- function(func, x.init, control, lower, upper) {
  ## By default, lower tolerance-- more likely to be met
  control <- modifyList(list(reltol=.Machine$double.eps^0.25,
                             parscale=rep(.1, length(x.init))),
                        control)

  check.bounds(lower, upper, x.init)
  if ( any(is.finite(lower) | is.finite(upper)) )
    func2 <- invert(boxconstrain(func, lower, upper))
  else
    func2 <- invert(func)

  ans <- subplex::subplex(x.init, func2, control)
  ans$value <- -ans$value
  names(ans)[names(ans) == "value"] <- "lnLik"

  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle (subplex): ",
            tolower(ans$message))
  
  ans
}

check.bounds <- function(lower, upper, x0=NULL) {
  if ( !is.null(x0) && (any(x0 < lower) || any(x0 > upper)) )
    stop("Starting parameter falls outside of problems bounds")
  if ( any(lower >= upper) )
    stop("'upper' must be strictly greater than 'lower'")
}

invert <- function(f) function(...) -f(...)

boxconstrain <- function(f, lower, upper, fail.value=-Inf) {
  function(x, ...) {
    if ( any(x < lower | x > upper) )
      fail.value
    else
      f(x, ...)
  }
}

set.defaults <- function(f, ..., defaults=NULL) {
  dots <- match.call(expand.dots=FALSE)[["..."]]
  if ( missing(defaults) )
    defaults <- dots
  else if ( is.list(defaults) )
    defaults <- c(dots, defaults)
  else
    stop("'defaults' must be a list")

  if ( is.null(defaults) )
    return(f)
  if ( !all(names(defaults) %in% names(formals(f))) )
    stop("Unknown defaults")
  att <- attributes(f)
  formals(f)[names(defaults)] <- defaults
  attributes(f) <- att[names(att) != "srcref"]
  f
}
