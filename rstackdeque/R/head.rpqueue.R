#' @export
#' @import utils
#' @title Return the head (front) of an rpqueue
#' 
#' @description Returns the first \eqn{n} elements of an rpqueue as an rpqueue, or all of the elements 
#' if \code{length(x) < n}.
#' 
#' @details Runs in \eqn{O(n)} time (in the size of the number of elements requested).
#' @param x rpqueue to get the head/top of.
#' @param n number of elements to get.
#' @param ... arguments to be passed to or from other methods (ignored).
#' @return an \code{\link{rpqueue}}.
#' @seealso \code{\link{rpqueue}}.
#' @examples 
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' q <- insert_back(q, "b")
#' q <- insert_back(q, "c")
#' 
#' qt <- head(q, n = 2)
#' print(qt)
head.rpqueue <- function(x, n = 6L, ...) {
  newrpqueue <- rpqueue()
  if(n < 0) {
    n = max(n, -1*length(x))
    n = length(x) + n
  } 
  if(n > length(x)) {
    n = length(x)
  } 
  if(n == 0) {
    return(newrpqueue)
  }
  for(i in seq(1,n)) {
    el <- peek_front(x)
    newrpqueue <- insert_back(newrpqueue, el)
    x <- without_front(x)
  }
  return(newrpqueue)
}