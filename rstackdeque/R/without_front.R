#' @export
#' @title Return a version of an rdeque or rpqueue without the front element
#' 
#' @details Simply returns a version of the given structure without the front element.
#' The original is left alone.
#' 
#' @details Runs in \eqn{O(1)}-amortized time for rdeque when used non-persistently (see documentation for 
#' \code{\link{rdeque}} for details), \eqn{O(1)} worst-case for rpqueue. Will
#' throw an error if the structure is empty to begin with.
#' 
#' @param x rdeque or rpqueue to remove elements from.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a version of the rdeque or rpqueue with the front element removed.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' d <- insert_back(d, "c")
#' 
#' d2 <- without_front(d)
#' print(d2)
#' 
#' d3 <- without_front(d2)
#' print(d3)
#' 
#' print(d)
without_front <- function(x, ...) {UseMethod("without_front", x)}