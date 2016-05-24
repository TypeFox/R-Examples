#' @export
#' @title Create an rstack pre-filled from a given input
#' 
#' @description Creates a new rstack from a given input. Coerces input to a 
#' list first using \code{as.list}, the element at \code{x[[1]]} becomes the top of the new rstack.
#'
#' @details Runs in \eqn{O(N)} in the size of the input. Because data frames return a list of 
#' columns when run through \code{as.list}, running \code{as.rstack} results in a stack of
#' columns, rather than a stack of rows.
#' @param x input to convert to a stack.
#' @param ... additional arguments to be passed to or from methods.
#' @return a new rstack.
#' @seealso \code{\link{rstack}}.
#' @examples
#' s <- as.rstack(1:20)
#' print(s)
#' 
#' s <- as.rstack(1:200000)
#' print(s)
#' 
#' ## A stack with only 5 elements, one for each column
#' oops <- as.rstack(iris)
#' print(oops)
as.rstack <- function(x, ...) {UseMethod("as.rstack", x)}