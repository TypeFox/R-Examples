#' @export
#' @title Return the number of elements in an rpqueue
#' 
#' @description Returns the number of elements in an rpqueue.
#' 
#' @details Runs in \eqn{O(1)} time, as this information is stored seperately and updated on every insert/remove.
#' @param x rpqueue to get the length of.
#' @return a vector of length 1 with the number of elements.
#' @seealso \code{\link{empty}} for checking whether an rpqueue is empty.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' print(length(q))         # 1
#' q <- insert_back(q, "b")
#' print(length(q))         # 2
#' @export
length.rpqueue <- function(x) {
  return(length(x$l) + length(x$r))
}