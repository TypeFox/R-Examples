#' @export
#' @title Assign to/modify the front of an rdeque or rpqueue
#' 
#' @description Allows modification access to the front of a deque or queue.
#' 
#' @details Runs in \eqn{O(1)} worst case time. Throws an error if the deque is empty.
#' @param x rdeque or rpqueue to modify the front element of.
#' @param value value to assign to the front data element.
#' @param ... additional arguments to be passed to or from methods.
#' @return modified rdeque or rpqueue.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' 
#' peek_front(d)$a <- 100
#' print(d)
#' 
#' peek_front(d) <- data.frame(a = 100, b = 100)
#' 
#' 
#' 
#' q <- rpqueue()
#' q <- insert_front(d, data.frame(a = 1, b = 1))
#' q <- insert_front(d, data.frame(a = 1, b = 1))
#' 
#' peek_front(q)$a <- 100
#' print(q)
#' 
#' peek_front(q) <- data.frame(a = 100, b = 100)
`peek_front<-` <- function(x, ..., value) { UseMethod("peek_front<-", x) }