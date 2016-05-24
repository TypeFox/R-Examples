#' @export
#' @title Assign to/modify the back of an rdeque
#' 
#' @description Allows modification access to the back of a deque.
#' 
#' @details Runs in \eqn{O(1)} worst case time. Throws an error if the deque is empty.
#' @param d rdeque to modify the back element of.
#' @param value value to assign to the back data element.
#' @param ... additional arguments to be passed to or from methods.
#' @return modified rdeque.
#' @seealso \code{\link{peek_back.rdeque}} for accessing the back element.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' 
#' peek_back(d)$a <- 100
#' print(d)
#' 
#' peek_back(d) <- data.frame(a = 100, b = 100)
`peek_back<-` <- function(d, ..., value) { UseMethod("peek_back<-", d) }