#' @export
#' @title Return the data element at the front of an rdeque
#' 
#' @description Simply returns the data element sitting at the front of the rdeque,
#' leaving the rdeque alone.
#' 
#' @details Runs in \code{O(1)} worst-case time.
#' @param x rdeque to peek at.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return data element existing at the front of the rdeque.
#' @seealso \code{\link{without_front}} for removing the front element.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' e <- peek_front(d)
#' print(e)
#' print(d)
#' 
#' ## Assigning to the front data element with peek_front:
#' d <- rdeque()
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' 
#' peek_front(d)$a <- 100
#' print(d)
#' 
#' peek_front(d) <- data.frame(a = 100, b = 100)
#' print(d)
peek_front.rdeque <- function(x, ...) {
  if(length(x) < 1) {
    stop("cannot peek_front() into a queue that is empty, try checking with empty() first")
  }
  if(length(x$l) > 0) {
    return(peek_top(x$l))
    # invariant: if l is empty but the deque is not, r has only one element
  } else {
    return(peek_top(x$r))
  }
}