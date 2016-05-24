#' @export
#' @title Fix an rdeque
#' @description A method used behind the scenes to provide \eqn{O(1)}-amortized time for most operations. 
#' Runs in \eqn{O(n)} time worst case; restructures the rdeque so that the two internal rstacks 
#' are roughly the same length.
#' @param d The rdeque to fix.
#' @param ... additional arguments to be passed to or from methods.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @return a fixed deque.
fixd.rdeque <- function(d, ...) {
  if(length(d) < 2 | (length(d$l) > 6 & length(d$r) > 6)) {
    return(d)
  } else {
    alllist <- c(as.list(d$l), rev(as.list(d$r)))
    mid <- as.integer(length(alllist)/2)
    left <- alllist[1:mid]
    right <- rev(alllist[(mid+1):length(alllist)])
    
    newd <- rdeque()
    newd$l <- as.rstack(left)
    newd$r <- as.rstack(right)
    return(newd)
  }
}