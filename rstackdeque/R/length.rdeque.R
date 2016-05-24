#' @export
#' @title Return the number of elements in an rdeque
#' 
#' @description Returns the number of elements in an rdeque.
#' 
#' @details Runs in \eqn{O(1)} time, as this information is stored seperately and updated on every insert/remove.
#' @param x rdeque to get the length of.
#' @return a vector of length 1 with the number of elements.
#' @seealso \code{\link{empty}} for checking whether an rdeque is empty.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' print(length(d))         # 1
#' d <- insert_back(d, "b")
#' print(length(d))         # 2
#' @export
length.rdeque <- function(x) {
  return(x$l$len + x$r$len)
}