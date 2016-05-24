#' Printing a Distributed Matrix
#' 
#' Print method for a distributed matrices.
#' 
#' Print method for class \code{ddmatrix}.
#' 
#' If argument \code{all=TRUE}, then a modified version of the ScaLAPACK TOOLS
#' routine PDLAPRNT is used to print the entire distributed matrix.  The matrix
#' will be printed in column-major fashion, with one element of the matrix per
#' line. If \code{all=FALSE} then the \code{name=} argument is ignored.
#' 
#' @param x,object
#' numeric distributed matrix
#' @param ... 
#' additional arguments
#' @param all 
#' control for whether the entire distributed matrix should be
#' printed to standard output
#' @param name 
#' character string that will be printed to standard output along
#' with the matrix elements
#' 
#' @return The function silently returns 0.
#' 
#' @examples
#' 
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:16, ncol=4)
#' dx <- as.ddmatrix(x) 
#' 
#' print(dx)
#' 
#' print(dx, all=T)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name ddmatrix-print
#' @rdname ddmatrix-print
NULL



#' @rdname ddmatrix-print
#' @export
setGeneric(name = "print", useAsDefault = base::print, package="pbdDMAT")



dmat.print <- function(dx)
{
  m <- dx@dim[1L]
  n <- dx@dim[2L]
  
  desca <- base.descinit(dim=dx@dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@ICTXT)
  
  base.rpdlaprnt(m=m, n=n, a=dx@Data, desca=desca)
}



#' @rdname ddmatrix-print
#' @export
setMethod("print", signature(x="ddmatrix"),
  function(x, ..., all=FALSE, name = "x"){
    if (all){
      assign(name, x)
      eval(parse(text = paste("dmat.print(", name, ")", sep = "") ))
    } else {
      ff <- paste(paste(format(base.firstfew(x, atmost=4), scientific=TRUE, digits=3), collapse=", "), ", ...", sep="")
      if (comm.rank()==0){
        grid <- base.blacs(x@ICTXT)
        cat(sprintf("\nDENSE DISTRIBUTED MATRIX\n---------------------------\n@Data:\t\t\t%s\nProcess grid:\t\t%dx%d\nGlobal dimension:\t%dx%d\n(max) Local dimension:\t%dx%d\nBlocking:\t\t%dx%d\nBLACS ICTXT:\t\t%d\n\n",
          ff, grid$NPROW, grid$NPCOL, x@dim[1], x@dim[2], x@ldim[1], x@ldim[2], x@bldim[1], x@bldim[2], x@ICTXT))
      }
    }
    
    pbdMPI::barrier()
    
    return( invisible(0) )
  }
)



#' @rdname ddmatrix-print
#' @export
setMethod("show", signature(object="ddmatrix"),
  function(object)
  {
    if (comm.rank()==0){
      grid <- base.blacs(object@ICTXT)
      cat(sprintf("\nDENSE DISTRIBUTED MATRIX\n---------------------------\nProcess grid:\t\t%dx%d\nGlobal dimension:\t%dx%d\n(max) Local dimension:\t%dx%d\nBlocking:\t\t%dx%d\nBLACS ICTXT:\t\t%d\n\n",
        grid$NPROW, grid$NPCOL, object@dim[1], object@dim[2], object@ldim[1], object@ldim[2], object@bldim[1], object@bldim[2], object@ICTXT))
    }
    
#    pbdMPI::barrier()
    
    return( invisible(0) )
  }
)

