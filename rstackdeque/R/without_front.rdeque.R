#' @export
#' @title Return a version of an rdeque without the front element
#' 
#' @description Simply returns a version of the given rdeque without the front element. Results in an error if the structure is empty.
#' The original rdeque is left alone.
#' 
#' @details Runs in \eqn{O(1)}-amortized time if the rdeque is used non-persistently (see documentation
#' of \code{\link{rdeque}} for details). If the given rdeque is empty, an error will be generated.
#' 
#' @param x rdeque to remove elements from.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return version of the rdeque with the front element removed.
#' @seealso \code{\link{insert_front}} for inserting elements.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' d <- insert_front(d, "c")
#' 
#' d2 <- without_front(d)
#' print(d2)
#' 
#' d3 <- without_front(d)
#' print(d3)
#' 
#' print(d)
without_front.rdeque <- function(x, ...) {
  ## if the length of l is 0, then r has only one element (invariant), so we can return an empty deque
  if(length(x$l) == 0) {
    return(rdeque())
  } 
  if(length(x) < 1) {
    stop("cannot run without_front() on an empty rdeque, check with empty() first")      
  }
  newd <- rdeque()
  newd$l <- x$l
  newd$r <- x$r
  newd$l <- without_top(newd$l)
  newd <- fixd(newd)
  return(newd)
}

