#' @export 
#' @title Return the data element at the front of an rdeque
#' 
#' @description Simply returns the data element sitting at the front of the deque,
#' leaving the deque alone.
#' 
#' @details Runs in \eqn{O(1)} worst-case time.
#' @param x rdeque to look at.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return data element at the front of the rdeque.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_back(d, "b")
#' e <- peek_front(d)
#' print(e)
#' print(d)
peek_front <- function(x, ...) {UseMethod("peek_front", x)}

