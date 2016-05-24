# This apply() operates on a MARGIN/ICTXT agreement, converting
# between data distributions as necessary.  This data movement
# is inarguably an inefficiency, but at present, it is the 
# best solution I can think of to preserve high level syntax
# on these data structures.

# Agreement occurs for ICTXT=1/MARGIN=2 and ICTXT=2/MARGIN=1

# For MARGIN=c(1,2), it is assumed that FUN is a one-to-one mapping 
# function which will not change the dimension of ddmatrix.



#' Apply Family of Functions
#' 
#' Apply a function to the margins of a distributed matrix.
#' 
#' 
#' If \code{reduce==TRUE} then a global matrix or vector (whichever is more
#' appropriate) will be returned. The argument \code{proc.dest=} behaves
#' exactly as in the \code{as.vector()} and \code{as.matrix()} functions of
#' \pkg{pbdDMAT}. If \code{reduce=FALSE} then a distributed matrix is returned.
#' Other acceptable arguments are \code{reduce="matrix"} and
#' \code{reduce="vector"} which demand global matrix or vector return,
#' respectively. This should generally be slightly more efficient than running
#' apply and then calling \code{as.vector()} or \code{as.matrix()}.
#' 
#' @param X 
#' distributed matrix
#' @param MARGIN 
#' subscript over which the function will be applied
#' @param FUN 
#' the function to be applied
#' @param ... 
#' additional arguments to FUN
#' @param reduce 
#' logical or string. See details
#' @param proc.dest 
#' Destination process (or 'all') if a reduction occurs
#' 
#' @return 
#' Returns a distributed matrix unless a reduction is requested, then a
#' global matrix/vector is returned.
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
#' y <- as.vector(apply(x, 1, mean))
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Extraction
#' @aliases apply
#' @name ddmatrix-apply
#' @rdname ddmatrix-apply
NULL

setGeneric(name = "apply", useAsDefault = base::apply, package="pbdDMAT")

#' @rdname ddmatrix-apply
#' @export
setMethod("apply", signature(X="ddmatrix"),
  function(X, MARGIN, FUN, ..., reduce=FALSE, proc.dest="all")
  {
    # idiot proofing
    if (missing(MARGIN))
      comm.stop('argument "MARGIN" is missing, with no default')
    else if (MARGIN != 1 && MARGIN != 2 && !all(MARGIN == c(1,2)))
      comm.stop('argument "MARGIN" must be 1, 2 or c(1,2) for a distributed matrix')
    
    oldCTXT <- X@ICTXT
    oldbldim <- X@bldim
    
    # Margin = c(1,2)
    if (all(MARGIN == c(1,2)))
    {
      resultOK <- FALSE
      if(ownany(X)){
        olddim <- dim(X@Data)
        X@Data <- base::apply(X=X@Data, MARGIN=MARGIN, FUN=FUN, ...)
        # check whether or not dimesion matches
        # currently it only supports one-to-one mapping function
        resultOK <- all(dim(X@Data) == olddim)
      }else{
        resultOK <- TRUE
      }
      resultOK <- pbdMPI::comm.all(resultOK)
      if(!resultOK) {
        comm.stop('apply only supports one-to-one mapping funcion when margin=c(1,2)')
      }
      return(X)
    }
    else if (MARGIN==1) # Row margin
    {
      if (X@ICTXT!=2)
      {
        fudge <- max(floor(X@bldim/base.blacs(X@ICTXT)$NPCOLS), 1)
        X <- dmat.reblock(dx=X, bldim=fudge, ICTXT=2)
      }
      
      if(ownany(X))
        tmp <- base::apply(X@Data, MARGIN=1, FUN=FUN, ...)
      else
        # it is unsafe to apply on X@Data if X does not hold any submatrix.
        tmp <- NULL
      
      ## The following block should be reconsidered since tmp can be NULL.
      if (!is.null(tmp) && is.list(tmp))
      {
        if (!all(sapply(tmp, is.numeric)))
          comm.stop("Error : list object contains non-numeric data")
        if (proc.dest=='all')
          return( pbdMPI::allgather(tmp) )
        else
          return( gather(tmp, rank.dest=proc.dest) )
      }

      else if (!is.null(tmp) && !is.matrix(tmp))
        dim(tmp) <- c(base::length(tmp), 1L)
      else if (!is.null(tmp) && is.matrix(tmp))
        # when apply on margin=1, the returned matrix should be transposed back
        tmp <- t(tmp) 
      
      # now we need to determine new ddmatrix global dimension.
      # row number remains same with old ddmatrix.
      # we need to get global column number which is the maxium of local new column number.
      if(ownany(X))
      {
        lcolnum <- dim(tmp)[2L]
      } else {
        # if X does not hold any submatrix, it is 0.
        lcolnum <- 0
      }
      gcolnum <- comm.max(lcolnum)
      
      if(ownany(X))
      {
        X@dim <- c(X@dim[1L], gcolnum)
        X@Data <- tmp
        X@ldim <- dim(X@Data)
      }
      else
      {
        X@dim <- c(X@dim[1L], gcolnum)
      }
      
    }
    else if (MARGIN==2) # Column margin
    {
      if (X@ICTXT!=1)
        fudge <- max(floor(X@bldim/base.blacs(X@ICTXT)$NPROWS), 1)
        X <- dmat.reblock(dx=X, bldim=X@bldim/2, ICTXT=1)
      
      if(ownany(X))
      {
        tmp <- base::apply(X@Data, MARGIN=2, FUN=FUN, ...)
      } 
      else 
      {
        # it is unsafe to apply on X@Data if X does not hold any submatrix.
        tmp <- NULL
      }
      
      ## The following block should be reconsidered since tmp can be NULL.
      if (!is.null(tmp) && is.list(tmp))
      {
        if (!all(sapply(tmp, is.numeric)))
          comm.stop("Error : list object contains non-numeric data")
        if (proc.dest=='all')
          return( pbdMPI::allgather(tmp) )
        else
          return( gather(tmp, rank.dest=proc.dest) )
      }
      else if (!is.null(tmp) && !is.matrix(tmp))
        dim(tmp) <- c(1L, base::length(tmp))
        
      # now we need to determine new ddmatrix global dimension.
      # column number remains same with old ddmatrix.
      # we need to get global row number which is the maxium of local row number of new submatrices.
      if(ownany(X))
        lrownum <- dim(tmp)[1L]
      else
        lrownum <- 0
      
      grownum <- comm.max(lrownum)
      
      if(ownany(X))
      {
        X@dim <- c(grownum, X@dim[2L])
        X@Data <- tmp
        X@ldim <- dim(X@Data)
      }
      else
        X@dim <- c(grownum, X@dim[2L])
    }
    
    if (reduce==TRUE)
    {
      if (MARGIN==1)
        if (X@dim[2L]==1)
          X <- as.vector(X, proc.dest=proc.dest)
        else
          X <- as.matrix(X, proc.dest=proc.dest)
        
      if (MARGIN==2)
        if (X@dim[1L]==1)
          X <- as.vector(X, proc.dest=proc.dest)
        else
          X <- as.matrix(X, proc.dest=proc.dest)
    }
    else if (reduce=="matrix")
      X <- as.matrix(X, proc.dest=proc.dest)
    else if (reduce=="vector")
      X <- as.vector(X, proc.dest=proc.dest)
    
    
    if (is.ddmatrix(X))
      if (X@ICTXT != oldCTXT)
        X <- dmat.reblock(dx=X, bldim=oldbldim, ICTXT=oldCTXT)
    
    if (is.ddmatrix(X) && length(MARGIN)==1 && MARGIN==1 && ncol(X)!=1) {
        # make the returned matrix has the same behaviour with R original apply function on margin = 1
        X <- t(X)
    }
    
    return(X)
  }
)
