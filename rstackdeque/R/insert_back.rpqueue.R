#' @export
#' @title Insert an element into the back of an rpqueue
#' 
#' @description Returns a version of the queue with the new element in the back position.
#' 
#' @details Runs in \eqn{O(1)} time worst-case. Does not modify the original. 
#' @param x rpqueue to insert onto.
#' @param e element to insert.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified version of the rpqueue.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' q <- insert_back(q, "b")
#' print(q)
#' 
#' q2 <- insert_back(q, "c")
#' print(q2)
#' print(q)
insert_back.rpqueue <- function(x, e, ...) {
  newq <- rpqueue()
  newq$l <- x$l
  newq$r <- insert_top(x$r, e)
  newq$lhat <- x$lhat
  return(makeequal(newq))
}

