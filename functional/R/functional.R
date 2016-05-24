##' Pre-specify a procedures named parameters, returning a new procedure.
##'
##' Thanks, Byron Ellis.
##' \url{https://stat.ethz.ch/pipermail/r-devel/2007-November/047318.html}
##' @param FUN the function to be curried
##' @param ... the determining parameters
##' @return A new function partially determined
##' @export
##' @examples
##' double <- Curry(`*`, e1=2)
##' stopifnot(double(4) == 8)
Curry <- function(FUN,...) {
  .orig = list(...);
  function(...) do.call(FUN,c(.orig,list(...)))
}

##' Lazy curry; thanks, Jamie!
##' <https://github.com/klutometis/R-functional/issues/1>
##' @inheritParams Curry
##' @export
##' @examples
##' # day is not defined; thanks, Jamie Folson.
##' CurryL(function(...) match.call(),
##'        x=5,
##'        y=as.Date(day))(z=as.Date(day,"%Y"))
CurryL <- function(FUN, ...){
  .curried <- as.list(match.call())[c(-1,-2)]
  function(...){
    .args <- as.list(match.call())[-1]
    eval(substitute(do.call(FUN,c(.curried,.args))))
  }}

##' Negate a function; borrowed from src/library/base/R/funprog.R for
##' pre-2.7 Rs.
##' @param f the function to be negated
##' @return The negated function
##' @examples
##' is.even <- function(a) a%%2 == 0
##' is.odd <- Negate(is.even)
##' stopifnot(Reduce(`&&`, Map(is.odd, c(1, 3, 5))))
Negate <- function(f)
  function(...) ! match.fun(f)(...)

##' Compose an arbitrary number of functions.
##'
##' My Happy Hacking keyboard gave out during the writing of this
##' procedure; moment of silence, please.
##' @param \dots the functions to be composed
##' @return A composed function
##' @export
##' @examples
##' car <- function(list) list[[1]]
##' cdr <- function(list) list[2:length(list)]
##' cadr <- Compose(cdr, car)
##' stopifnot(cadr(c(1,2,3)) == 2)
Compose <- function(...) {
  fs <- list(...)
  
  ## Thanks, Matthew Lungberg.
  if (!all(sapply(fs, is.function)))
    stop("Argument is not a function")

  function(...) Reduce(function(x, f) f(x),
                       fs,
                       ...)
}

##' Identity function.
##'
##' Is concatenation benign?
##' @param \dots tautological arguments
##' @return The tautologized arguments, concatenated
##' @export
##' @examples
##' list.copy <- function(list)
##'   Reduce(Identity, list)
##' 
##' list <- c(1, 2, 3)
##' stopifnot(list.copy(list) == list)
Identity <- function(...) c(...)

##' Thanks, Gabor; see <http://stackoverflow.com/a/23726989>: swaps
##' the first two arguments in a function.
##' @param f The function whose arguments to swap
##' @return A function with swapped arguments
Swap <- function(f) function(x, y, ...) f(y, x, ...)

##' Composition with multiple arguments.
##'
##' Thanks, Alexander Davis!
##' @param \dots the functions to be composed
##' @return A composed function
##' @export
##' @examples
##' f <- function(x, y) x+y
##' g <- function(x) x*2
##' stopifnot(multi.argument.Compose(f, g)(1,1) == 4)
multi.argument.Compose <- function (...) {
    fs <- list(...)
    if (!all(sapply(fs, is.function)))
        stop("Argument is not a function")
    function(...) Reduce(function(x, f) f(x), fs[-1], fs[[1]](...))
}
