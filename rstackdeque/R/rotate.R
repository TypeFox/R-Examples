#' @export
#' @title Generic maintenance function for rpqueues
#' @description Generic maintenance function for rpqueues, called automatically when needed by other functions.
#' @details See \emph{Simple and Efficient Purely Functional Queues and Deques}, 
#' Okasaki 1995, J. Functional Programming, 5(4) 583 to 592 for information.
#' @param rpqueue rpqueue to rotate.
#' @param acclazylist lazy list accumulator.
#' @param ... additional arguments to be passed to or from methods (ignored).
#' @return a fully rotated rpqueue, but with the l list delayed in evaluation.
#' @seealso \code{\link{makeequal}} helper function that calls this one.
#' @references Okasaki, Chris. Purely Functional Data Structures. Cambridge University Press, 1999.
#' @description The lazy list accumulator, keeping the queue partially rotated. 
rotate <- function(rpqueue, acclazylist, ...) {UseMethod("rotate", rpqueue)}

