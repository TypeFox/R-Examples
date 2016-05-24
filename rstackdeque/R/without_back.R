#' @export
#' @title Return a version of an rdeque without the back element
#' 
#' @description Simply returns a version of the given rdeque without the back element
#' The original rdeque is left alone.
#' 
#' @details Runs in \eqn{O(1)}-amortized time if the rdeque is used non-persistently (see documentation
#' of \code{\link{rdeque}} for details). If the given rdeque is empty, an error will be generated.
#' 
#' @param d rdeque to remove elements from.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return version of the rdeque with the back element removed.
#' @seealso \code{\link{insert_back}} for inserting elements.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' d <- insert_front(d, "c")
#' 
#' d2 <- without_back(d)
#' print(d2)
#' 
#' d3 <- without_back(d)
#' print(d3)
#' 
#' print(d)
without_back <- function(d, ...) {UseMethod("without_back", d)}