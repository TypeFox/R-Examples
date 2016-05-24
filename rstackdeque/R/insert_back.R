#' @export
#' @title Insert an element into the back of an rdeque or rpqueue
#' 
#' @description Returns a version of the deque/queue with the new element in the back position.
#' 
#' @details Runs in \eqn{O(1)} time worst-case. Does not modify the original. 
#' @param x rdeque or rpqueue to insert onto.
#' @param e element to insert.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified version of the rdeque or rpqueue.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' d <- rdeque()
#' d <- insert_back(d, "a")
#' d <- insert_back(d, "b")
#' print(d)
#' 
#' d2 <- insert_back(d, "c")
#' print(d2)
#' print(d)
insert_back <- function(x, e, ...) {UseMethod("insert_back", x)}