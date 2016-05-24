#' Non-Distributed object to Distributed Object Converters
#' 
#' A simplified interface to the \code{distribute()} and \code{redistribute()}
#' functions.
#' 
#' A simplified wrapper for the \code{distribute()} function, especially in the
#' case that the matrix \code{x} is global (which you really should not ever
#' let happen outside of testing, but I won't stop you).
#' 
#' The function will only work if \code{x} is stored on all processes, or
#' \code{x} is stored on a single process (does not matter which) and every
#' other process has NULL stored for x.
#' 
#' If several processes own pieces of the matrix \code{x}, then you can not use
#' this function. You will have to create an appropriate \code{ddmatrix} on all
#' processes and redistriubte the data with the \code{redistribute()} function.
#' 
#' As usual, the \code{ICTXT} number is the BLACS context corresponding to the
#' process grid onto which the output distributed matrix will be distributed.
#' 
#' @param x 
#' a numeric matrix
#' @param ...
#' Additional arguments.
#' @param bldim 
#' the blocking dimension for block-cyclically distributing the
#' matrix across the process grid.
#' @param xCTXT 
#' the BLACS context number for initial distribution of the matrix
#' \code{x}.
#' @param ICTXT 
#' BLACS context number for return.
#' 
#' @return Returns a distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' if (comm.rank()==0){
#'   x <- matrix(1:16, ncol=4)
#' } else {
#'   x <- NULL
#' }
#' 
#' dx <- as.ddmatrix(x, bldim=2)
#' dx
#' 
#' ### Can also be common to all ranks
#' y <- matrix(1:25, 5, bldim=2)
#' dy <- as.ddmatrix(y)
#' dy
#' 
#' finalize()
#' }
#' 
#' @keywords Distributing Data
#' @name as.ddmatrix
#' @rdname as.ddmatrix
setGeneric(name="as.ddmatrix", 
  function(x, ...) 
    standardGeneric("as.ddmatrix"), 
  package="pbdDMAT"
)



base.mat.to.ddmat <- function(x, bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
{
  if (!is.matrix(x))
    comm.stop("input 'x' must be a matrix") 
  else if (length(bldim) == 1) 
    bldim <- rep(bldim, 2) 
  else if (diff(bldim) != 0)
    comm.warning("Most ScaLAPACK routines do not allow for non-square blocking.  This is highly non-advised.")
  
  dim <- dim(x)
  ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
  descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
  
  out <- base.mksubmat(x=x, descx=descx)
  
  dx <- new("ddmatrix", Data=out, dim=dim, ldim=dim(out), bldim=bldim, ICTXT=ICTXT)
  
  return(dx)
}



# distribute a matrix from process (0,0) to the full ICTXT grid
#' @rdname as.ddmatrix
#' @export
distribute <- function(x, bldim=.pbd_env$BLDIM, xCTXT=0, ICTXT=.pbd_env$ICTXT)
{
  if (length(bldim)==1)
    bldim <- rep(bldim, 2L)
  
  if (!is.matrix(x) && is.null(x)){
    x <- matrix(0)
    iown <- FALSE
  } else
    iown <- TRUE
  
  if (iown)
    dim <- dim(x)
  else
    dim <- c(0, 0)
  
  ldim <- dim(x)
  
  if (!is.double(x))
    storage.mode(x) <- "double"
  
  blacs_ <- blacs(xCTXT)
  
  if (blacs_$NPROW > 1)
    dim[1] <- pbdMPI::allreduce(dim[1], op='sum')
  else
    dim[1] <- pbdMPI::allreduce(dim[1], op='max')
  
  if (blacs_$NPCOL > 1)
    dim[2] <- pbdMPI::allreduce(dim[2], op='sum')
  else
    dim[2] <- pbdMPI::allreduce(dim[2], op='max')
  
  if (all(ldim==0))
    ldim <- c(1,1)
  
  dx <- new("ddmatrix", Data=x, dim=dim, ldim=ldim, bldim=dim, ICTXT=xCTXT)
  
  if (xCTXT != ICTXT)
    dx <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=ICTXT)
  else if (any(dx@bldim != bldim))
    dx <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=dx@ICTXT)
  
  return( dx )
}



# Distribute dense, in-core matrices
dmat.as.ddmatrix <- function(x, bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
{
  nprocs <- pbdMPI::comm.size()
  owns <- pbdMPI::allreduce(is.matrix(x), op='sum')
  
  # owned by one process 
  if (owns==1)
  { 
    iown <- is.matrix(x)
    if (iown)
      iown <- pbdMPI::comm.rank()
    else
      iown <- 0
    iown <- pbdMPI::allreduce(iown, op='max')
    return( distribute(x=x, bldim=bldim, xCTXT=0, ICTXT=ICTXT) )
  } 
  # global ownership is assumed --- this should only ever really happen in testing
  else if (owns==nprocs)
    return( base.mat.to.ddmat(x, bldim=bldim, ICTXT=ICTXT) )
  # neither of these two cases
  else 
    comm.stop("Matrix 'x' is defined on some, but not all processes. Consider using the redistribute() function.")
}



#' @rdname as.ddmatrix
#' @export
setMethod("as.ddmatrix", signature(x="matrix"), 
  dmat.as.ddmatrix
)



#' @rdname as.ddmatrix
#' @export
setMethod("as.ddmatrix", signature(x="NULL"), 
  dmat.as.ddmatrix
)



#' @rdname as.ddmatrix
#' @export
setMethod("as.ddmatrix", signature(x="vector"), 
  function(x, bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
    dmat.as.ddmatrix(matrix(x), bldim=bldim, ICTXT=ICTXT)
)
