#' @export
#' @title Create a new empty rpqueue
#' 
#' @description Creates a new, empty, rpqueue ready for use.
#' 
#' @return a new rpqueue.
#' @details An rpqueue provides "First In, First Out" (FIFO) access; envisaged
#' as a queue, elements may be inserted at the back and removed from the front. Unlike
#' \code{\link{rdeque}}, access is gauranteed \eqn{O(1)} worst case even when used 
#' persistently, though in most situations rdeques will be faster in practice
#' (see the documentation for \code{\link{rdeque}} for details). 
#' 
#' Other handy functions
#' include \code{as.list} and \code{as.data.frame} (the latter of which requires that
#' all elements can be appended to become rows of a data frame in a reasonable manner). 
#' @seealso \code{\link{rstack}} for a fast LIFO-only structure.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @examples
#' q <- rpqueue()
#' q <- insert_back(q, "a")
#' q <- insert_back(q, "b")
#' print(q)
#' 
#' q2 <- without_front(q)
#' print(q2)
#' print(q)
#' 
#' b <- peek_front(q)
#' print(b)
rpqueue <- function() {
  newq <- new.env(parent = emptyenv())
  newq$lhat <- rstack()
  newq$l <- rstack()
  newq$r <- rstack()
  class(newq) <- "rpqueue"
  return(newq)
}