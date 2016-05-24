#' @export
#' @title Return the data element at the top of an rstack
#' 
#' @description Simply returns the data element sitting at the top of the rstack,
#' leaving the rstack alone.
#' 
#' @details Runs in \code{O(1)} worst-case time.
#' @param s rstack to peek at.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return data element existing at the top of the rstack.
#' @seealso \code{\link{without_top}} for removing the top element.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' e <- peek_top(s)
#' print(e)
#' print(s)
#' 
#' ## Assigning to the top data element with peek_top:
#' s <- rstack()
#' s <- insert_top(s, data.frame(a = 1, b = 1))
#' s <- insert_top(s, data.frame(a = 1, b = 1))
#' 
#' peek_top(s)$a <- 100
#' print(s)
#' 
#' peek_top(s) <- data.frame(a = 100, b = 100)
peek_top.rstack <- function(s, ...) {
  if(is.null(s$head)) {
    stop("cannot peek() at the top of an empty stack, try checking with empty() first")
  } else {
    return(s$head$data)
  }
}