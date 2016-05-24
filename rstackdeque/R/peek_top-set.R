#' @export
#' @title Assign to/modify the top of an rstack
#' 
#' @description Allows modification access to the top of a stack.
#' 
#' @details Runs in \eqn{O(1)} worst case time. Throws an error if the rstack is \code{\link{empty}}. Modifies the element in place (i.e., is not side-effect-free).
#' @param s rstack to modify the first element of.
#' @param value value to assign to the top data element.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return modified rstack.
#' @seealso \code{\link{peek_top}} for accessing the top data element.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, data.frame(a = 1, b = 1))
#' s <- insert_top(s, data.frame(a = 1, b = 1))
#' 
#' peek_top(s)$a <- 100
#' print(s)
#' 
#' peek_top(s) <- data.frame(a = 100, b = 100)
`peek_top<-` <- function(s, ..., value) { UseMethod("peek_top<-", s) }