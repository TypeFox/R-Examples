#' eigen
#' 
#' Eigenvalue decomposition for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' @param symmetric 
#' logical, if \code{TRUE} then the matrix is assumed to be
#' symmetric and only the lower triangle is used.  Otherwise \code{x} is
#' inspected for symmetry.
#' @param only.values 
#' logical, if \code{TRUE} then only the eigenvalues are
#' returned.  Otherwise both eigenvalues and eigenvectors are returned.
#' @param EISPACK
#' Ignored.
#' 
#' @return 
#' \code{eigen()} computes the eigenvalues, and eigenvectors if requested.  As
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- ddmatrix(1:9, 3, bldim=2)
#' 
#' y <- eigen(crossprod(x))
#' y
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Linear Algebra
#' @aliases eigen
#' @name ddmatrix-eigen
#' @rdname ddmatrix-eigen
NULL

setGeneric(name = "eigen", useAsDefault = base::eigen, package="pbdDMAT")

#' @rdname ddmatrix-eigen
#' @export
setMethod("eigen", signature(x="ddmatrix"), 
  function(x, symmetric, only.values=FALSE, EISPACK=FALSE)
  {
    if (x@dim[1L] != x@dim[2L])
      comm.stop("non-square matrix in 'eigen'")
    
    if (missing(symmetric)) 
      symmetric <- isSymmetric(x)
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    if (symmetric){
      if (only.values)
        jobz <- 'N'
      else
        jobz <- 'V'
      
      out <- base.rpdsyev(jobz=jobz, uplo='L', n=x@dim[2L], a=x@Data, desca=desca, descz=desca)
    } else {
      if (!only.values) # FIXME
        comm.stop("Currently only possible to recover eigenvalues from a non-symmetric matrix")
      
#      out <- base.pdgeeig(dx@Data, descx=desca)
    }
    
    return( out )
  }
)

