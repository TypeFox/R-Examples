# ------------------------------------------------
# PDTRAN:  Matrix transpose
# ------------------------------------------------

#' rpdtran
#' 
#' Transpose.
#' 
#' For advanced users only.
#' 
#' @param a
#' Matrix.
#' @param desca,descc
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdtran <- function(a, desca, descc)
{
  m <- descc[3L]
  n <- descc[4L]
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  
  cldim <- base.numroc(descc[3L:4L], descc[5L:6L], ICTXT=descc[2L])
  
  ret <- .Call(R_PDTRAN,
                as.integer(m), as.integer(n),
                a, as.integer(desca),
                as.integer(cldim), as.integer(descc))
  
  return(ret)
}

# ------------------------------------------------
# PDGEMM:  Matrix-Matrix multiplication
# ------------------------------------------------

#' rpdgemm
#' 
#' Matrix-Matrix Multiply.
#' 
#' For advanced users only.
#' 
#' @param transx,transy
#' 'T' or 'N' for transpose or not.
#' @param x,y
#' Matrix.
#' @param descx,descy,descc
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgemm <- function(transx, transy, x, descx, y, descy, descc)
{
  transx <- toupper(transx)
  transy <- toupper(transy)
  
  m <- descx[3L]
  n <- descy[4L]
  k <- descy[3L]
  
  cldim <- base.numroc(descc[3:4], descc[5:6], ICTXT=descc[2])
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  if (!is.double(y))
    storage.mode(y) <- "double"
  
  ret <- .Call(R_PDGEMM,
                transx, transy,
                as.integer(m), as.integer(n), as.integer(k),
                x, as.integer(descx),
                y, as.integer(descy),
                as.integer(cldim), as.integer(descc))
  
  return( ret )
}

