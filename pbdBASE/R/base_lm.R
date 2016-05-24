#' rpdgels
#' 
#' Linear model fitter via rank-revealing QR (with pivoting).
#' 
#' For advanced users only.
#' 
#' @param tol 
#' Numerical tolerance for the QR.
#' @param m,n
#' Problem size.
#' @param nrhs
#' Number of right hand sides.
#' @param a
#' Left hand side.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param b
#' Right hand side.
#' @param descb
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgels <- function(tol, m, n, nrhs, a, desca, b, descb)
{
#  # FIXME adjustment for weird lda issue
#  mxlda <- pbdMPI::allreduce(desca[9], op='max')
#  mxldb <- pbdMPI::allreduce(descb[9], op='max')
#  
#  if (desca[9]==1)
#    desca[9] <- mxlda
#  if (descb[9]==1)
#    descb[9] <- mxldb
  
  if (!is.double(a))
    storage.mode(a) <- "double"
  if (!is.double(b))
    storage.mode(b) <- "double"
  
  ltau <- as.integer(min(m, n))
  
  ret <- .Call(R_PDGELS,
            TOL=as.double(tol), M=as.integer(m), N=as.integer(n), NRHS=as.integer(nrhs),
            A=a, DESCA=as.integer(desca), B=b, DESCB=as.integer(descb),
            LTAU=ltau)
  
  
  # Sometimes R mistakenly frees these matrices...
  if (!base.ownany(dim=desca[5L:6L], bldim=desca[5L:6L], ICTXT=desca[2L]))
    ret$A <- matrix(0.0, 1L, 1L)
  
#  if (!base.ownany(dim=descb[5L:6L], bldim=descb[5L:6L], ICTXT=descb[2L]))
#  {
#    ret$EFF <- ret$RSD <- ret$FT <- matrix(0.0, 1L, 1L)
#  }
  
  
  if (ret$INFO!=0)
    comm.warning(paste("ScaLAPACK returned INFO=", ret$INFO, "; returned solution is likely invalid", sep=""))
  
  return( ret )
}



