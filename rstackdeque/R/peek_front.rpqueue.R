#' @export
#' @title Return the data element at the front of an rpqueue
#' 
#' @description Simply returns the data element sitting at the front of the rpqueue,
#' leaving the queue alone.
#' 
#' @details Runs in \code{O(1)} worst-case time.
#' @param x rpqueue to peek at.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return data element existing at the front of the queue.
#' @seealso \code{\link{without_front}} for removing the front element.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' q <- insert_back(q, "b")
#' e <- peek_front(q)
#' print(e)
#' print(q)
#' 
#' ## Assigning to the front data element with peek_front:
#' q <- rpqueue()
#' q <- insert_back(q, data.frame(a = 1, b = 1))
#' q <- insert_back(q, data.frame(a = 1, b = 1))
#' 
#' peek_front(q)$a <- 100
#' print(q)
#' 
#' peek_front(q) <- data.frame(a = 100, b = 100)
#' print(q)
peek_front.rpqueue <- function(x, ...) {
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