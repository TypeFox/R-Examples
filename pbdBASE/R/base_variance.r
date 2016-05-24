#' Column Variances
#' 
#' Computes the variances of a ScaLAPCK-like distributed matrix.
#' Significantly faster than using \code{apply()}, even in compared
#' to the performance differences you would find comparing these
#' two approaches using just base R.
#' 
#' @param x
#' The matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' 
#' @export
base.pdclvar <- function(x, descx)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_PDCLVAR, x, as.integer(descx), as.integer(dim(x)[2L]))
  
  return( ret )
}
