#' @export
#' @title Convert an rpqueue to a list
#' 
#' @description Converts an rpqueue to a list, where the front of the queue becomes
#' the first element of the list, the second-from-front the second, and so on. 
#' 
#' @details Runs in \eqn{O(N)} time in the size of the rpqueue, but the generated list is pre-allocated for efficiency.
#' @param x rpqueue to convert.
#' @param ... additional arguments passed to as.list after initial conversion to list.
#' @return a list containing the elements of the rpqueue in front-to-back order.
#' @seealso \code{\link{as.data.frame.rpqueue}} and the generic \code{\link{as.list}}.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' q <- insert_back(q, "b")
#' 
#' qlist <- as.list(q)
#' print(qlist)
#' @export
as.list.rpqueue <- function(x, ...) {
  return(as.list(c(as.list(x$l), rev(as.list(x$r))), ...))
}

