#' @export
#' @import utils
#' @title Return the first n elements of an rdeque as an rdeque
#' 
#' @description Returns the first n elements of a deque as a deque, or all of the elements if its length is less than n.
#' 
#' @details Runs in \eqn{O(n)} time (in the size of the number of elements requested).
#' @param x rdeque to get the head of.
#' @param ... arguments to be passed to or from other methods (ignored).
#' @param n number of elements to get.
#' @return a new rdeque.
#' @examples 
#' d <- rdeque()
#' d <- insert_back(d, "a")
#' d <- insert_back(d, "b")
#' d <- insert_back(d, "c")
#' 
#' dt <- head(d, n = 2)
#' print(dt)
head.rdeque <- function(x, n = 6L, ...) {
  newdeque <- rdeque()
  if(n < 0) {
    n = max(n, -1*length(x))
    n = length(x) + n
  } 
  if(n > length(x)) {
    n = length(x)
  } 
  if(n == 0) {
    return(newdeque)
  }
  for(i in seq(1,n)) {
    newdeque <- insert_back(newdeque, peek_front(x))
    x <- without_front(x)
  }
  return(newdeque)
}