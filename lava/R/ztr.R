##' Trace operator
##'
##' Calculates the trace of a square matrix.
##' @param A Square numeric matrix
##' @return \code{numeric}
##' @author Klaus K. Holst
##' @seealso \code{\link{crossprod}}, \code{\link{tcrossprod}}
##' @keywords math algebra
##' @examples
##'
##' tr(diag(1:5))
##'
##' @export
`tr` <-
function(A) {
  if (length(A)==1)
    return(A)
  if(!is.matrix(A))
    stop("argument of 'tr' should be a matrix.")
  n <- nrow(A)
  if (!n)
    stop("0 x 0 matrix")
  if (n != ncol(A))
    stop("non-square matrix")
  if (any(!is.finite(A)))
    stop("infinite or missing values")
  return(sum(diag(A)))
}
