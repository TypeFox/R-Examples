#' @export
#' @title Return the number of elements in an rstack
#' 
#' @description Returns the number of elements in an rstack.
#' 
#' @details Runs in \eqn{O(1)} time, as this information is stored seperately and updated on every insert/remove.
#' @param x rstack to get the length of.
#' @return a vector of length 1, which the number of elements of the stack.
#' @seealso \code{\link{empty}} for checking whether an rstack is empty.
#' @examples
#' s <- rstack()
#' s <- insert_top(s, "a")
#' print(length(s))         # 1
#' s <- insert_top(s, "b")
#' print(length(s))         # 2
length.rstack <- function(x) {
  return(x$len)
}