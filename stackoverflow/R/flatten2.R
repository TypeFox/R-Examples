#' Flatten a list without type coercion
#'
#' @section Changed Feb 19, 2015 by njf:
#' Rather than calculating length, preallocate more than needed. 
#'
#' @param x a nested list
#' @param len guess of output length
#'
#' @author \href{http://stackoverflow.com/users/662787/tommy}{Tommy}, 
#'   \href{http://stackoverflow.com/users/271616/joshua-ulrich}{Joshua Ulrich}, 
#'   \href{http://stackoverflow.com/users/980833/josh-obrien}{Josh O'Brien}, 
#'   Neal Fultz
#' @references \url{http://stackoverflow.com/questions/8139677/how-to-flatten-a-list-to-a-list-without-coercion}
#' @export

flatten2 <- function(x, len=1024) {
  y <- vector('list', len)
  i <- 1L
  rapply(x, function(x) { y[[i]] <<- x; i <<- i + 1L })
  y[seq_len(i - 1L)]
}