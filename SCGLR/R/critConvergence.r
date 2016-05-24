#' @title Auxiliary function for controlling SCGLR fitting
#' @description Auxiliary function for \code{scglr} fitting used to 
#' construct a convergence control argument.
#' @export 
#' @param tol positive convergence threshold.
#' @param maxit integer, maximum number of iterations.
#' @return a list containing elements named as the arguments.
critConvergence <- function (tol = 1e-6, maxit = 50) {
  if (!is.numeric(tol) || tol <= 0) 
    stop("value of 'tolerance' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  list(tol = tol, maxit = maxit)
}
