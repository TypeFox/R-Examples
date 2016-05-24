#' @export
#' @title Insert an element into the back of an rdeque
#' 
#' @description Returns a version of the deque with the new element in the back position.
#' 
#' @details Runs in \eqn{O(1)} time worst-case. Does not modify the original. 
#' @param x rdeque to insert onto.
#' @param e element to insert.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified version of the rdeque.
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
insert_back.rdeque <- function(x, e, ...) {
  newd <- rdeque()
  newd$r <- insert_top(x$r, e)
  newd$l <- x$l
  newd <- fixd(newd)
  return(newd)
}

