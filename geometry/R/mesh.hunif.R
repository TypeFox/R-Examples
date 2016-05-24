##' Uniform desired edge length function of position to allow easy
##' definition of regions when passed as the \code{fh} argument of
##' \code{\link{distmesh2d}} or \code{\link{distmeshnd}}.
##'
##' @title Uniform desired edge length
##' @param p A \code{n}-by-\code{m} matrix, each row representing a
##' point in an \code{m}-dimensional space.
##' @param ... additional arguments (not used)
##' @return Vector of ones of length \code{n}.
##' @author Raoul Grasman; translated from original Matlab sources of Per-Olof
##' Persson.
##' @seealso \code{\link{distmesh2d}} and \code{\link{distmeshnd}}.
##' @export
mesh.hunif <- function (p, ...) {
  if (!is.matrix(p)) 
    stop("Input `p' should be matrix.")
  return(rep(1, nrow(p)))
}
