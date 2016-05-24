#' crossprod
#' 
#' Crossproduct.
#' 
#' For advanced users only.
#' 
#' @param uplo
#' Triangle whose values to use.
#' @param trans
#' tcrossprod or crossprod.
#' @param x
#' Matrix to crossprod.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param descc
#' ScaLAPACK descriptor array of output.
#' 
#' @export
base.crossprod <- function(uplo, trans, x, descx, descc)
{
  trans <- toupper(trans)
  uplo <- toupper(uplo)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call(R_PDCROSSPROD,
                  uplo, trans, x, as.integer(descx),
                  as.integer(cldim), as.integer(descc))
  
  return( ret )
}



#' pdchtri
#' 
#' Inverse of cholesky.
#' 
#' For advanced users only.
#' 
#' @param uplo
#' Triangle whose values to use.
#' @param x
#' Matrix to crossprod.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param descc
#' ScaLAPACK descriptor array of output.
#' 
#' @export
base.pdchtri <- function(uplo, x, descx, descc)
{
  # FIXME move row/col adjustment down to Fortran (currently in calling DMAT chol2inv method)
  
  uplo <- toupper(uplo)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  ret <- .Call(R_PDCHTRI, 
                uplo, x, as.integer(dim(x)), as.integer(descx), 
                as.integer(cldim), as.integer(descc))
  
  return( ret )
}

