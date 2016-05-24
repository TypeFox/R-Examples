#' @export
#' @title Create a new empty rdeque
#' 
#' @description Creates a new, empty, rdeque ready for use.
#' 
#' @return a new empty rdeque.
#' @details An rdeque provided both "Last In, First Out" (LIFO) and "First In, First Out" (FIFO) access;
#' envisaged as queue, elements may be added or removed from the front or the back with \code{\link{insert_front}},
#' \code{\link{insert_back}}, \code{\link{without_front}}, and \code{\link{without_back}}. The front and back 
#' elements may be retrieved with \code{\link{peek_front}} and \code{\link{peek_back}}.
#' 
#' Internally, rdeques are stored as a pair of rstacks; they provide \eqn{O(1)}-amortized insertion and removal,
#' so long as they are not used persistently (that is, the variable storing the deque is always replaced
#' by the modified version, e.g. \code{s <- without_front(s)}). When used persistently (e.g. \code{s2 <- without_front(s)}, which results
#' in two deques being accessible), this cannot be gauranteed. When an \eqn{O(1)} worst-cast, fully
#' persistent FIFO queue is needed, the rpqueue from this package provides these semantics. However,
#' there is a constant-time overhead for rpqueues, such that rdeques are often more efficient (and flexible)
#' in practice, even in when used persistently.
#' 
#' Other useful functions
#' include \code{as.list} and \code{as.data.frame} (the latter of which requires that
#' all elements can be appended to become rows of a data frame in a reasonable manner).
#' @seealso \code{\link{rstack}} for a fast LIFO-only structure.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' d <- rdeque()
#' d <- insert_front(d, "a")
#' d <- insert_front(d, "b")
#' d <- insert_back(d, "c")
#' d <- insert_back(d, "d")
#' print(d)
#' 
#' d2 <- without_back(d)
#' print(d2)
#' print(d)
#' 
#' b <- peek_front(d)
#' print(b)
#' @export
rdeque <- function() {
  d <- new.env(parent = emptyenv())
  d$l <- rstack()
  d$r <- rstack()
  class(d) <- "rdeque"
  return(d)
}