## function_utils.R
##   - Utility functions for R functions
##
## RGP - a GP system for R
## 2010 Oliver Flasch (oliver.flasch@fh-koeln.de)
## with contributions of Thomas Bartz-Beielstein, Olaf Mersmann and Joerg Stork
## released under the GPL v2
##

##' Create a new function stub
##'
##' Creates and returns a new function stub without capturing any environment variables.
##'
##' @note Always use this function to dynamically generate new functions that are not clojures
##' to prevent hard to find memory leaks.
##' @param envir The new function closure's environment, defaults to \code{globalenv()}.
##' @return A new function that does not take any arguments and always returns \code{NULL}.
new.function <- function(envir = globalenv()) {
  # The following line has to be inside this function to prevent capturing the lexical
  # environment, which would cause a hard to find memory leak:
  fun <- function() NULL
  # If you use the last line outside of this function, by shure to initialize the
  # environment of the generated function object to the global environment, like done here:
  environment(fun) <- envir
  fun 
}

##' Create a new R closure given a function body expression and an argument list 
##'
##' Creates a R closure (i.e. a function object) from a body expression and an argument
##' list. The closure's environment will be the default environment.
##'
##' @param fbody The function body, given as a R expression.
##' @param fargs The formal arguments, given as a list or vector of strings.
##' @param envir The new function closure's environment, defaults to \code{globalenv()}.
##' @return A formal argument list, ready to be passed via \code{\link{formals}}.
##' @export
makeClosure <- function(fbody, fargs, envir = globalenv())
  .Call("make_closure", fbody, fargs, envir)

##' Create a new function argument list from a list or vector of strings
##'
##' Creates a formal argument list from a list or vector of strings, ready to be assigned via
##' \code{\link{formals}}.
##'
##' @param fargs The formal arguments, given as a list or vector of strings.
##' @return A formal argument list, ready to be passed via \code{\link{formals}}.
new.alist <- function(fargs) {
  alistargs <- Reduce(function(a,b) paste(a,b,sep="=,") , fargs, "", right = TRUE)
  alistargslen <- nchar(alistargs)
  if (alistargslen != 0)
    alistargs <- substr(alistargs, 1, alistargslen-1)
  eval(parse(text = paste("alist(", alistargs, ")", sep="")), baseenv())
}

##' A variant of do.call that ignores unused arguments
##'
##' @param what What to call (either a function or a character vector naming a function in
##'   \code{envir}.
##' @param args The args for the call, these may include arguments not used by \code{what}.
##' @param quote Whether to quote the arguments.
##' @param envir The environment within which to evaluate the call.
##' @return The result of the call.
do.call.ignore.unused.arguments <- function(what, args, quote = FALSE, envir = parent.frame()) {
  fcn <- if (is.function(what)) {
    what
  } else if (is.character(what) && exists(what, mode = "function", envir = envir)) {
    get(what, mode = "function", envir = envir)
  } else stop("do.call.ignore.unused.arguments: The parameter 'what' must be a function or a character string naming a function.")
  fcnFormalNames <- names(formals(fcn))
  okFcnFormalNames <- intersect(names(args), fcnFormalNames)
  okFcnFormals <- args[okFcnFormalNames]
  do.call(fcn, okFcnFormals, quote = quote, envir = envir)
}

##' Repeatedly apply a function
##'
##' Repeatedly apply a function \code{f} to an argument \code{arg}, additional arguments \code{...}
##' are supplied unchanged in each call. E.g. \code{iterate(3, foo, 42.14, "bar")} is equivalent to
##' \code{foo(foo(foo(42.14, "bar"), "bar"), "bar")}.
##'
##' @param n The number of times to apply \code{f}, must be >= 0. If 0, \code{arg} is returned.
##' @param f The function to apply.
##' @param arg The argument to repeatedly apply \code{f} to.
##' @param ... Additional argument to pass to \code{f} at each application.
##' @return The result of repeatedly applying \code{f}.
iterate <- function(n, f, arg, ...) {
  stopifnot(n >= 0)
  result <- arg
  if (n > 0) for (i in 1:n) result <- f(result, ...) # TODO is this the most efficient way to iterate ?
  result
}

##' Determine the number of arguments of a function
##'
##' Tries to determine the number of arguments of function.
##'
##' @param f The function to determine the arity for.
##' @return The arity of the function \code{f}.
arity <- function(f) {
  if (is.primitive(f)) {
    arity.primitive(f)
  } else if (is.function(f)) {
    length(formals(f))
  } else {
    stop("could not determine arity of object \"", f, "\" (of class ", class(f), ")")
  }
}

##' Determine the number of arguments of a primitive function
##'
##' Tries to determine the number of arguments of a primitive R function by lookup in a
##' builtin table.
##'
##' @param f The primitive to determine the arity for.
##' @return The arity of the primitive \code{f}.
arity.primitive <- function(f) {
  if (identical(f, `+`) ||
      identical(f, `*`) ||
      identical(f, `-`) ||
      identical(f, `/`) ||
      identical(f, `^`) ||
      identical(f, `%%`) ||
      identical(f, `%/%`) ||
      identical(f, `<`) ||
      identical(f, `<=`) ||
      identical(f, `==`) ||
      identical(f, `!=`) ||
      identical(f, `>=`) ||
      identical(f, `>`) ||
      identical(f, `&`) ||
      identical(f, `|`) ||
      identical(f, `&&`) ||
      identical(f, `||`) ||
      identical(f, xor))
    2
  else if (identical(f, log) ||
           identical(f, logb) ||
           identical(f, atan2))
    2
  else if (identical(f, log10) ||
           identical(f, log2) ||
           identical(f, log1p) ||
           identical(f, exp) ||
           identical(f, expm1) ||
           identical(f, sin) ||
           identical(f, cos) ||
           identical(f, tan) ||
           identical(f, asin) ||
           identical(f, acos) ||
           identical(f, atan) ||
           identical(f, abs) ||
           identical(f, sqrt) ||
           identical(f, `!`))
    1
  else
    stop("could not determine arity of primitive")
}

##' Tabulate an n-ary function
##'
##' Creates a data frame of values for the n-ary function \code{f} at the sample
##' locations given in \code{...}.
##'
##' @param f The function to tabulate.
##' @param ... For each dimension, a vector of sample points to calculate \code{f} at.
##' @return A data frame of function values of \code{f}.
##' @export 
tabulateFunction <- function(f, ...) {
  xs <- expand.grid(...)
  colnames(xs) <- NULL
  ys <- apply(xs, 1, function(args) do.call(f, as.list(args)))
  functionTable <- cbind(xs, ys)
  colnames(functionTable) <- c(paste("x", 1:length(list(...)), sep = ""), "y")
  functionTable
} 

