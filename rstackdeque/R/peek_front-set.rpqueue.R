#' @export
#' @title Assign to/modify the front of an rpqueue
#' 
#' @description Allows modification access to the front of a queue.
#' 
#' @details Runs in \eqn{O(1)} worst case time. Throws an error if the rpqueue is \code{\link{empty}}. Modifies the element in place (i.e., is not side-effect-free).
#' @param x rpqueue to modify the front element of.
#' @param value value to assign to the front data element.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified rpqueue.
#' @seealso \code{\link{peek_front.rpqueue}} for accessing the front data element.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, data.frame(a = 1, b = 1))
#' q <- insert_back(q, data.frame(a = 1, b = 1))
#' 
#' peek_front(q)$a <- 100
#' print(q)
#' 
#' peek_front(q) <- data.frame(a = 100, b = 100)
#' print(q)
`peek_front<-.rpqueue` <- function(x, ..., value) {
  if(length(x) < 1) {
    stop("cannot assign to the front of an empty deque or queue, try checking with empty() first")
  }
  
  if(length(x$l) > 0) {
    peek_top(x$l, ...) <- value 
    # invariant: if l is empty but the deque is not, r has only one element
  } else {
    peek_top(x$r, ...) <- value
  }
  
  return(x)
}

