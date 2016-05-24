#' @export
#' @title Check if an rstack, rdeque, or rpqueue is empty
#' 
#' @description Check if an rstack, rdeque, or rpqueue is empty.
#' @details Runs in \eqn{O(1)} time for all types.
#' 
#' @param x rstack, rdeque, or rpqueue to check.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return logical vector of length 1.
#' @examples
#' s <- rstack()
#' print(empty(s))        ## TRUE
#' s <- insert_top(s, "a")
#' print(empty(s))        ## FALSE
empty <- function(x, ...) {UseMethod("empty", x)}