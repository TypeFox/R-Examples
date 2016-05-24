#' @export
#' @title Convert an rdeque to a list
#' 
#' @description Converts an rdeque to a list, where the front of the deque becomes
#' the first element of the list, the second-from-front the second, and so on. 
#' 
#' @details Runs in \eqn{O(N)} time in the size of the rdeque, but the generated list is pre-allocated for efficiency.
#' @param x rdeque to convert.
#' @param ... additional arguments passed to as.list after initial conversion to list.
#' @return a list containing the elements of the rdeqeue in front-to-back order.
#' @seealso \code{\link{as.data.frame.rstack}} and the generic \code{\link{as.list}}.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' 
#' dlist <- as.list(d)
#' print(dlist)
as.list.rdeque <- function(x, ...) {
  return(as.list(c(as.list(x$l), rev(as.list(x$r))), ...))
}