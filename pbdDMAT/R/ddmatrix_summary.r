#' Distributed Matrix Summary
#' 
#' Summarize a distributed matrix.  Gives min, max, mean, etc. by column.
#' 
#' The return is on process 0 only.
#' 
#' @param object 
#' numeric distributed matrix
#' @param ...
#' Additional arguments.
#' 
#' @return 
#' A table on processor 0, \code{NULL} on all other processors.
#' 
#' @examples
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
#' summary(dx)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name ddmatrix-summary
#' @rdname ddmatrix-summary
NULL



#' @rdname ddmatrix-summary
#' @export
setGeneric(name = "summary", useAsDefault = base::summary, package="pbdDMAT")



#' @rdname ddmatrix-summary
#' @export
setMethod("summary", signature(object="ddmatrix"),
  function(object)
  {
    if (object@ICTXT != 1)
    {
      newbldim <- c(object@dim[1], ceiling(object@bldim[2] / object@dim[2]))
      object <- dmat.redistribute(object, bldim=newbldim, ICTXT=1)
    }
    
    if (ownany(object))
    {
      lret <- summary(object@Data)
      
      ret_names <- sapply(X=1:object@ldim[2], FUN=
        function(i) 
          base.l2g_coord(ind=c(1, i), dim=object@dim, bldim=object@bldim, ICTXT=1)[2]
      )
    } 
    else 
    {
      lret <- NULL
      ret_names <- NULL
    }
    
    ret <- gather(lret)
    ret_names <- gather(ret_names)
    
    if (comm.rank() == 0)
    {
      ret <- ret[which(!sapply(ret, is.null))]
      ret <- array(unlist(ret), c(6L, object@dim[2L]))
      row.names(ret) <- rep("", 6L)
      
      ret_names <- ret_names[which(!sapply(ret_names, is.null))]
      ret_names <- unlist(ret_names)
      
      print(ret_names)
      
      colnames(ret) <- paste("V", ret_names, sep="")
      
      if (any(colnames(ret) != ret_names))
        ret <- ret[, paste("V", 1L:object@dim[2L], sep=""), drop=F]
      
      return( as.table(ret) )
    }
    else
      return( invisible(NULL) )
  }
)
