#' p_matpow_by_squaring_wrap
#' 
#' Matrix power by squaring.
#' 
#' For advanced users only.
#' 
#' @param A
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param b
#' Power.
#' 
#' @export
base.p_matpow_by_squaring_wrap <- function(A, desca, b=1)
{
  b <- as.integer(b)
  desca <- as.integer(desca)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call(R_p_matpow_by_squaring, A, desca, b)
  
  return( ret )
}



#' p_matexp_pade_wrap
#' 
#' Pade' expansion.
#' 
#' For advanced users only.
#' 
#' @param A
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param p
#' Order of the Pade' approximation.
#' 
#' @export
base.p_matexp_pade_wrap <- function(A, desca, p=6)
{
  desca <- as.integer(desca)
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  ret <- .Call(R_p_matexp_pade, A, desca, as.integer(p))
  
  return( ret )
}

