#' @export
#' @title Assign to/modify the front of an rdeque
#' 
#' @description Allows modification access to the front of a deque.
#' 
#' @details Runs in \eqn{O(1)} worst case time. Throws an error if the rdeque is \code{\link{empty}}. Modifies the element in place (i.e., is not side-effect-free).
#' @param x rdeque to modify the front element of.
#' @param value value to assign to the front data element.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified rdeque.
#' @seealso \code{\link{peek_front.rdeque}} for accessing the front data element.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' 
#' peek_front(d)$a <- 100
#' print(d)
#' 
#' peek_front(d) <- data.frame(a = 100, b = 100)
#' print(d)
`peek_front<-.rdeque` <- function(x, ..., value) {
  if(length(x) < 1) {
    stop("cannot assign to the front of an empty deque, try checking with empty() first")
  }
  x$l$head$data <- value
  return(x)
}