#' Distributed object to Matrix Converters
#' 
#' Converts a distributed matrix into a non-distributed matrix.
#' 
#' The \code{proc.dest=} argument accepts either the BLACS grid position or the
#' MPI rank if the user desires a single process to own the matrix.
#' Alternatively, passing the default value of \code{'all'} will result in all
#' processes owning the matrix. If only a single process owns the undistributed
#' matrix, then all other processes store \code{NULL} for that object.
#' 
#' @param x 
#' numeric distributed matrix
#' @param ...
#' Additional arguments.
#' @param proc.dest 
#' destination process for storing the matrix
#' @param attributes 
#' logical, specifies whether or not the current attributes
#' should be preserved.
#' 
#' @return 
#' Returns an ordinary R matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' dx <- ddmatrix(1:16, ncol=4, bldim=2)
#' y <- as.matrix(dx, proc.dest=0)
#' 
#' comm.print(y)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods
#' @name as.matrix
#' @rdname as.matrix
NULL


setGeneric(name = "as.matrix", useAsDefault = base::as.matrix, package="pbdDMAT")



# create a global matrix from a ddmatrix
dmat.gmat <- function(dx, proc.dest="all")
{
  xattrs <- attributes(dx@Data)
  names <- xattrs$dimnames
  
  ICTXT <- dx@ICTXT
  
  dim <- dx@dim
  ldim <- dx@ldim
  bldim <- dx@bldim
  
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  if (any(dim==0)){
    if (proc.dest[1L] == "all" || proc.dest==comm.rank())
      out <- matrix(nrow=dim[1], ncol=dim[2])
    else
      out <- NULL
    return(out)
  }
  
  if (proc.dest[1]=='all')
    rsrc <- csrc <- -1
  else {
#    dest <- base.pcoord(ICTXT=ICTXT, PNUM=proc.dest)
#    rsrc <- dest[[1]]
#    csrc <- dest[[2]]
    rsrc <- proc.dest[1]
    csrc <- proc.dest[2]
  }
  
  out <- base.mkgblmat(dx@Data, descx=descx, rsrc=rsrc, csrc=csrc)
  
  if (is.null(out))
    return(out)
  else {
    if (length(xattrs)>1){
      if (length(names)>0)
        xattrs$dimnames <- NULL
      xattrs$dim <- dim(out)
      attributes(out) <- xattrs
    }
    
    return( out )
  }
}



# Undistribute a distributed matrix --- ONLY to be used in testing
base.as.matrix <- function(x, proc.dest="all") 
{
  if (proc.dest=='all'){
    ret <- dmat.gmat(dx=x, proc.dest="all")
    return( ret )
  }
  else if (is.numeric(proc.dest)){
    if (base::length(proc.dest)==1){
      blacs_ <- base.blacs(x@ICTXT)
      if (pbdMPI::comm.rank()==proc.dest)
        proc.dest <- c(blacs_$MYROW, blacs_$MYCOL)
      else
        proc.dest <- c(0, 0)
      proc.dest <- pbdMPI::allreduce(proc.dest, op='max')
    } 
    else if (base::length(proc.dest)>2)
      comm.stop("Invalid destination process 'proc.dest'")
    
    ret <- dmat.gmat(dx=x, proc.dest=proc.dest)
    return( ret )
  }
  
  comm.stop("Invalid destinaction process 'proc.dest'")
}



#' @rdname as.matrix
#' @export
setMethod("as.matrix", signature(x="ddmatrix"), 
  function(x, proc.dest="all", attributes=TRUE)
  {
    # convert ddmatrix attributes too
    if (attributes){
      ddms <- sapply(attributes(x@Data), is.ddmatrix)
      if (any(ddms)){
        for (att in which(ddms)){
          if (any(attributes(x@Data)[[att]]@ldim == 1)){
            attributes(x@Data)[[att]] <- as.vector(attributes(x@Data)[[att]])
          }
          else
            attributes(x@Data)[[att]] <- as.matrix(attributes(x@Data)[[att]])
        }
      }
    }
    
    
    ret <- base.as.matrix(x=x, proc.dest=proc.dest)
    
    if (is.logical(x@Data))
      storage.mode(ret) <- "logical"
    
    return( ret )
  }
)



##' @rdname as.matrix
##' @export
#setMethod("as.matrix", signature(x="dmat"),
#  function(x)
#  {
#    mat <- matrix(0.0, x@dim[1L], x@dim[2L])
#    
#    dim <- x@dim
#    nrows <- dim[1L]
#    
#    nrows.local <- dmat_ldim(nrows)
#    ldim <- c(nrows.local, dim[2L])
#    
#    start <- dmat_index(nrows)
#    end <- start + nrows.local - 1L
#    
#    if (ldim[1L] > 0)
#      mat[start:end, ] <- x@Data
#    
#    # FIXME make this bcast later, too lazy atm
#    mat <- allreduce(mat)
#    
#    return( mat )
#  }
#)



##' @rdname as.matrix
##' @export
#setMethod("as.matrix", signature(x="dsmatrix"),
#  function(x)
#  {
#    y <- as.matrix(as.dmat(x))
#    
#    return( y )
#  }
#)


