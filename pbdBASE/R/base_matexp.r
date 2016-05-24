#' matexp
#' 
#' Serial matrix exponentiation.
#' 
#' For advanced users only.
#' 
#' @param A
#' Matrix to exponentiate.
#' @param p
#' Pade' expansion size.
#' @param t
#' Scaling factor.
#' 
#' @export
base.matexp <- function(A, p=6, t=1)
{
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  R <- .Call(R_matexp, A, as.integer(p), as.double(t))
  
  return( R )
}

