#' (Un)Distribute
#' 
#' (Un)Distribute matrix.
#' 
#' For advanced users only.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' 
#' @rdname lclgblmat
#' @export
base.mksubmat <- function(x, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L], fixme=TRUE)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  subx <- .Call(R_MKSUBMAT, x, as.integer(ldim), as.integer(descx))
  
  return( subx )
}



#' @param rsrc,csrc
#' Row/column source.
#' @rdname lclgblmat
#' @export
base.mkgblmat <- function(x, descx, rsrc, csrc)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_MKGBLMAT, 
       x, as.integer(descx), as.integer(rsrc), as.integer(csrc))
  
  return( ret )
  
}



#' dallreduce
#' 
#' Allreduce
#' 
#' For advanced users only.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param op
#' Operation.
#' @param scope
#' Rows, columns, or both.
#' 
#' @export
base.dallreduce <- function(x, descx, op='sum', scope='All')
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_DALLREDUCE, 
        x, as.integer(dim(x)), as.integer(descx), as.character(op), as.character(scope))
  
  return( ret )
}



#' tri2zero
#' 
#' Zero Triangle
#' 
#' For advanced users only.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param uplo
#' Triangle.
#' @param diag
#' Zero diagonal as well.
#' 
#' @export
base.tri2zero <- function(x, descx, uplo='L', diag='N')
{
  uplo <- toupper(uplo)
  diag <- toupper(diag)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  ret <- .Call(R_PTRI2ZERO, 
               uplo, diag, x, as.integer(dim(x)), as.integer(descx))
  
  return( ret )
}



#' pdsweep
#' 
#' Matrix-Vector Sweep
#' 
#' For advanced users only.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param vec
#' Vector
#' @param MARGIN
#' Rows or columns.
#' @param FUN
#' Function.
#' 
#' @export
base.pdsweep <- function(x, descx, vec, MARGIN, FUN)
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (!is.double(vec))
    storage.mode(vec) <- "double"
  
  ret <- .Call(R_PDSWEEP, 
               x, as.integer(dim(x)), as.integer(descx), vec, as.integer(length(vec)), as.integer(MARGIN), as.character(FUN))
  
  return( ret )
}



#' diag
#' 
#' Grab diagonal or create distributed diagonal matrix.
#' 
#' For advanced users only.
#' 
#' @param x
#' Matrix.
#' @param descx
#' ScaLAPACK descriptor array.
#' @param proc.dest
#' Who owns the result.
#' 
#' @name diag
#' @rdname diag
#' @export
base.ddiagtk <- function(x, descx, proc.dest='all')
{
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  if (proc.dest[1L] == 'all')
    rdest <- cdest <- -1
  else {
    if (length(proc.dest)==1){
      src <- base.pcoord(ICTXT=descx[2L], PNUM=proc.dest)
      rsrc <- src[[1L]]
      csrc <- src[[2L]]
    }
  }
  
  ldiag <- min(descx[3L:4L])
  
  ret <- .Call(R_PDGDGTK, 
               x, as.integer(dim(x)), as.integer(descx), as.integer(ldiag),
               as.integer(rdest), as.integer(cdest))
  
  return( ret )
}

#' @param diag
#' Diagonal.
#' 
#' @rdname diag
#' @export
base.ddiagmk <- function(diag, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L])
  
  if (!is.double(diag))
    storage.mode(diag) <- "double"
  
  out <- .Call(R_PDDIAGMK, 
               as.integer(ldim), as.integer(descx), diag, as.integer(length(diag)))
  
  return( out )
}



#' dhilbmk
#' 
#' Create Hilbert matrix.
#' 
#' For advanced users only.
#' 
#' @param n
#' Size.
#' 
#' @export
base.dhilbmk <- function(n)
{
  n <- as.integer(n)
  
  ret <- .Call(R_DHILBMK, n)
  
  return( ret )
}



#' pdhilbmk
#' 
#' Create Hilbert matrix.
#' 
#' For advanced users only.
#' 
#' @param descx
#' ScaLAPACK descriptor matrix.
#' 
#' @export
base.pdhilbmk <- function(descx)
{
  descx <- as.integer(descx)
  ldim <- as.integer(base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L], fixme=TRUE))
  
  ret <- .Call(R_PDHILBMK, ldim, descx)
  
  return( ret )
}



#' pdmkcpn1
#' 
#' Create Companion Matrix
#' 
#' For advanced users only.
#' 
#' @param coef
#' Coefficients vector.
#' @param descx
#' ScaLAPACK descriptor array.
#' 
#' @export
base.pdmkcpn1 <- function(coef, descx)
{
  ldim <- base.numroc(dim=descx[3L:4L], bldim=descx[5L:6L], ICTXT=descx[2L])
  
  if (!is.double(coef))
    storage.mode(coef) <- "double"
  
  out <- .Call(R_PDMKCPN1, as.integer(ldim), as.integer(descx), coef)
  
  return( out )
}



