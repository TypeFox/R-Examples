#' @export
#' @title Return the data element at the back of an rdeque
#' 
#' @description Simply returns the data element sitting at the back of the rdeque,
#' leaving the rdeque alone.
#' 
#' @details Runs in \code{O(1)} worst-case time.
#' @param d rdeque to peek at.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return data element existing at the back of the rdeque.
#' @seealso \code{\link{without_back}} for removing the front element.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' e <- peek_back(d)
#' print(e)
#' print(d)
#' 
#' ## Assigning to the front data element with peek_front:
#' d <- rdeque()
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' d <- insert_front(d, data.frame(a = 1, b = 1))
#' 
#' peek_back(d)$a <- 100
#' print(d)
#' 
#' peek_back(d) <- data.frame(a = 100, b = 100)
#' print(d)
peek_back <- function(d, ...) {UseMethod("peek_back", d)}