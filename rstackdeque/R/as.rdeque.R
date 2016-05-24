#' @export
#' @title Create a pre-filled rdeque from a given input
#' 
#' @description Creates a new rdeque from a given input. Coerces input to a 
#' list first using \code{as.list}, the element at \code{x[[1]]} becomes the front of the new rdeque.
#'
#' @details Runs in \eqn{O(N)} in the size of the input. Because data.frames return a list of 
#' columns when run through \code{as.list}, running \code{as.rdeque} results in a deque of
#' columns, rather than a deque of rows.
#' @param x input to convert to an rdeque.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a new rdeque.
#' @seealso \code{\link{rdeque}}.
#' @examples
#' d <- as.rdeque(1:20)
#' print(d)
#' 
#' d <- as.rdeque(1:200000)
#' print(d)
#' 
#' ## A deck with only 5 elements, one for each column
#' oops <- as.rdeque(iris)
#' print(oops)
as.rdeque <- function(x, ...) {UseMethod("as.rdeque", x)}