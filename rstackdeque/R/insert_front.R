#' @export
#' @title Insert an element into the front of an rdeque
#' 
#' @description Returns a version of the deque with the new element in the front position.
#' 
#' @details Runs in \eqn{O(1)} time worst-case. Does not modify the original rdeque. 
#' @param d rdeque to insert onto.
#' @param e element to insert.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified version of the rdeque.
#' @seealso \code{\link{without_front}} for removing the front element.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' print(d)
#' 
#' d2 <- insert_front(d, "c")
#' print(d2)
#' print(d)
insert_front <- function(d, e, ...) {UseMethod("insert_front", d)}