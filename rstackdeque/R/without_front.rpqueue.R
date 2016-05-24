#' @export
#' @title Return a version of an rpqueue without the front element
#' 
#' @description Simply returns a version of the given rpqueue without the front element. Results in an error if the structure is empty.
#' The original rpqueue is left alone.
#' 
#' @details Runs in \eqn{O(1)} worst case time. 
#' 
#' @param x rpqueue to remove elements from.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return version of the rpqueue with the front element removed.
#' @seealso \code{\link{peek_front}} for accessing the front element.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' q <- insert_back(q, "b")
#' q <- insert_back(q, "c")
#' 
#' q2 <- without_front(q)
#' print(q2)
#' 
#' q3 <- without_front(q)
#' print(q3)
#' 
#' print(q)
without_front.rpqueue <- function(x, ...) {
  if(length(x) < 1) {
    stop("cannot run without_front() on an empty rpqueue, try checking with empty() first")
  }
  newq <- rpqueue()
  newq$l <- without_top(x$l)
  newq$r <- x$r
  newq$lhat <- x$lhat
  return(makeequal(newq))
}