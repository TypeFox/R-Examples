#' Extract or construct a diagonal matrix.
#'
#' This works like \code{\link{diag}} except when \code{x} is a single
#' integer value.  If \code{x} is a single integer value then it
#' assumes that you want a 1 by 1 matrix with the value set to \code{x}
#' @param x a matrix, vector or 1D array, or missing.
#' @param nrow,ncol optional dimensions for the result when \code{x} is not a matrix.
#' @return matrix with diagonal elements set to \code{x}
#' @export
#' @seealso \code{\link{diag}}
#' @keywords array
#' @examples
#' diag(5)
#' diagn(5)

diagn <- function(x=1, nrow=length(x), ncol=nrow)
  if(is.matrix(x)) diag(x=x) else diag(x=x, nrow=nrow, ncol=ncol)
