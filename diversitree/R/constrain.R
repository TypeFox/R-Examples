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
