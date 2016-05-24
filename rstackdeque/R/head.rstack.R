#' @export
#' @import utils
#' @title Return the head (top) of an rstack
#' 
#' @description Returns the top \eqn{n} elements of an rstack as an stack, or all of the elements 
#' if \code{length(x) < n}.
#' 
#' @details Runs in \eqn{O(n)} time (in the size of the number of elements requested).
#' @param x rstack to get the head/top of.
#' @param n number of elements to get.
#' @param ... arguments to be passed to or from other methods (ignored).
#' @return an \code{\link{rstack}}.
#' @seealso \code{\link{rstack}}.
#' @examples 
#' s <- rstack()
#' s <- insert_top(s, "a")
#' s <- insert_top(s, "b")
#' s <- insert_top(s, "c")
#' 
#' st <- head(s, n = 2)
#' print(st)
#' print(s)
head.rstack <- function(x, n = 6L, ...) {
  newstack <- rstack()
  if(n < 0) {
    n = max(n, -1*length(x))
    n = length(x) + n
  } 
  if(n > length(x)) {
    n = length(x)
  } 
  if(n == 0 | n < -1*length(x)) {
    return(newstack)
  }
  node <- x$head
  for(i in seq(1,n)) {
    if(!is.null(node)) {
      newstack <- insert_top(newstack, node$data)
      node <- node$nextnode
    }
  }
  return(rev(newstack))
}

