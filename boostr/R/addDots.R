#' @title Extend a function's signature to include '...'
#'
#' @description modifies the input function to accept extra arguments, 
#' if it doesn't already.
#'
#' @param func the function whose signature is to be modified.
#' @param .verbose logical flag indicating whether warnings should be displayed
#' or not.
#' 
#' @return a function with the same body as \code{func} but whose signature now
#' includes \code{...}.
#' 
#' @examples
#' \dontrun{
#' f <- function(x,y) x^2 + y^2
#' 
#' g <- addDots(f)
#' g
#' 
#' h <- addDots(g, .verbose=TRUE)
#' }
addDots <- function(func, .verbose=FALSE) {
  sigOld <- formals(func)
  if ("..." %in% names(sigOld)) {
    if (.verbose) warning("'...' is already in function signature!")
    return(func)
  } else {
    sigNew <- c(sigOld, alist("..." = ))
    formals(func) <- sigNew
    return(func)
  }
}