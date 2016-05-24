#' @export
#' @title Generic maintenance function for rpqueues
#' @description Generic maintenance function for rpqueues, called automatically when needed by other functions.
#' @details See \emph{Simple and Efficient Purely Functional Queues and Deques}, 
#' Okasaki 1995, J. Functional Programming, 5(4) 583 to 592 for information.
#' @param rpqueue rpqueue to makeequal.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a "fixed" rpqueue.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @seealso \code{\link{rotate}} helper function that calls this one.
makeequal <- function(rpqueue, ...) {UseMethod("makeequal", rpqueue)}

