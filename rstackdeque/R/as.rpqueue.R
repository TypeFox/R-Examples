#' @export
#' @title Create a pre-filled rpqueue from a given input
#' 
#' @description Creates a new rpqueue from a given input. Coerces input to a 
#' list first using \code{as.list}, the element at \code{x[[1]]} becomes the front of the new queue.
#'
#' @details Runs in \eqn{O(N)} in the size of the input. Because data.frames return a list of 
#' columns when run through \code{as.list}, running \code{as.rpqueue} results in a queue of
#' columns, rather than a queue of rows.
#' @param x input to convert to an rpqueue.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a new rpqueue.
#' @seealso \code{\link{rpqueue}}.
#' @examples
#' d <- as.rpqueue(1:20)
#' print(d)
#' 
#' ## A queue with only 5 elements, one for each column
#' oops <- as.rdeque(iris)
#' print(oops)
as.rpqueue <- function(x, ...) {UseMethod("as.rpqueue", x)}