### General redistribution
# redistribute data from one context and/or blocking factor to another
dmat.reblock <- function(dx, bldim=dx@bldim, ICTXT=.pbd_env$ICTXT)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  if (length(bldim)==1)
    bldim <- rep(bldim, 2)
  
  dim <- dx@dim
  m <- dim[1]
  n <- dim[2]
  xattrs <- attributes(dx@Data)
  
  ldimB <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
  TldimB <- ldimB # true ldimB
  
  # lda's of 1 infuriate pdgemr2d
  mxx <- pbdMPI::allreduce(max(dx@ldim), op='max')
  mxb <- pbdMPI::allreduce(max(ldimB), op='max')
  
  if (all(dx@ldim==1))
    dx@ldim[1] <- mxx
  if (all(ldimB==1))
    ldimB[1] <- mxb
  
  if (!ownany(x=dx))
    ldimB <- c(0,0)
  
#  if (dx@dim[1]>1 && pbdMPI::allreduce(dx@ldim[1], op='max')==1)
#    dx@ldim[1] <- mxx
#  if (pbdMPI::allreduce(ldimB[1], op='max')==1)
#    ldimB[1] <- mxb
  
  descx <- base.descinit(dim=dim, bldim=dx@bldim, ldim=dx@ldim, ICTXT=dx@ICTXT)
  descy <- base.descinit(dim=dim, bldim=bldim, ldim=ldimB, ICTXT=ICTXT)
  
  if (!is.double(dx@Data))
    storage.mode(dx@Data) <- "double"
  
  ret <- base.rpdgemr2d(x=dx@Data, descx=descx, descy=descy)
  
  dy <- new("ddmatrix", Data=ret, dim=dim, ldim=TldimB, bldim=bldim, ICTXT=ICTXT)
  
  if (length(xattrs) > 1){
    xattrs$dim <- dy@ldim
    attributes(dy@Data) <- xattrs
  }
  
  
  return( dy )
}

### TODO fix this mess
dmat.redistribute <- dmat.reblock



#' Distribute/Redistribute matrices across the process grid
#' 
#' Takes either an R matrix and distributes it as a distributed matrix, or
#' takes a distributed matrix and redistributes it across a (possibly) new
#' BLACS context, using a (possibly) new blocking dimension.
#' 
#' \code{distribute()} takes an R matrix \code{x} stored on the processes in
#' some fashion and distributes it across the process grid belonging to
#' \code{ICTXT}. If a process is to call \code{distribute()} and does not yet
#' have any ownership of the matrix \code{x}, then that process should store
#' \code{NULL} for \code{x}.
#' 
#' How one might typically use this is to read in a non-distributed matrix on
#' the first process, store that result as the R matrix \code{x}, and then have
#' the other processes store \code{NULL} for \code{x}. Then calling
#' \code{distribute()} returns the distributed matrix which was distributed
#' according to the options \code{bldim} and \code{ICTXT}.
#' 
#' Using an \code{ICTXT} value other than zero is not recommended unless you
#' have a good reason to. Use of other such contexts should only be considered
#' for advanced users, preferably those with knowledge of ScaLAPACK.
#' 
#' \code{redistribute()} takes a distributed matrix and redistributes it to the
#' (possibly) new process grid with BLACS context \code{ICTXT} and with the
#' (possibly) new blocking dimension \code{bldim}. The original BLACS context
#' is \code{dx@@ICTXT} and the original blocking dimension is \code{dx@@bldim}.
#' 
#' These two functions are essentially simple wrappers for the ScaLAPACK
#' function PDGEMR2D, with the above described behavior. Of note, for
#' \code{distribute()}, \code{dx@@ICTXT} and \code{ICTXT} must share at least one
#' process in common. Likewise for \code{redistribute()} with \code{xCTXT} and
#' \code{ICTXT}.
#' 
#' Very general redistributions can be done with \code{redistribute()}, but
#' thinking in these terms is an acquired skill.  For this reason, several
#' simple interfaces to this function have been written.
#' 
#' @param dx 
#' numeric distributed matrix
#' @param bldim 
#' the blocking dimension for block-cyclically distributing the
#' matrix across the process grid.
#' @param ICTXT 
#' BLACS context number for return.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @keywords BLACS Distributing Data
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
#' dx <- distribute(x, bldim=c(4,4))
#' print(dx)
#' 
#' dx <- redistribute(dx, bldim=c(3,3))
#' print(dx)
#' 
#' 
#' finalize()
#' }
#' 
#' @name redistribute
#' @rdname redistribute
#' @export
reblock <- dmat.reblock

#' @rdname redistribute
#' @export
redistribute <- dmat.reblock




#' Distribute/Redistribute matrices across the process grid
#' 
#' Takes either an R matrix and distributes it as a distributed matrix, or
#' takes a distributed matrix and redistributes it across a (possibly) new
#' BLACS context, using a (possibly) new blocking dimension.
#' 
#' These functions are simple wrappers of the very general
#' \code{redistribute()} funciton.  Different
#' distributed matrix distributions of note can be classified into three
#' categories: block, cyclic, and block-cyclic.
#' 
#' \code{as.block()} will convert \code{ddmatrix} into one which is merely
#' "block" distributed, i.e., the blocking factor is chosen in such a way that
#' there will be no cycling.  By default, this new blocking factor will be
#' square.  This can result in some raggedness (some processors owning less
#' than others --- or nothing) if the matrix is far from square itself.
#' However, the methods of factoring \code{ddmatrix} objects, and therefore
#' anything that relies on (distributed) matrix factorizations such as
#' computing an inverse, least squares solution, etc., require that blocking
#' factors be square.  The matrix will not change BLACS contexts.
#' 
#' \code{as.rowblock()} will convert a distributed matrix into one which is
#' distributed by row into a block distributed matrix.  That is, the rows are
#' stored contiguously, and different processors will own different rows, but
#' with no cycling.  In other words, it block redistributes the data across
#' context 2.
#' 
#' \code{as.colblock()} is the column-wise analogue of \code{as.rowblock()}.
#' In other words, it block redistributes the data across context 1.
#' 
#' \code{as.rowcyclic()} is a slightly more general version of
#' \code{as.rowblock()}, in that the data will be distributed row-wise, but
#' with the possibility of cycling, as determined by the blocking factor.  In
#' other words it block-cyclically redistributes the data across context 2.
#' 
#' \code{as.colcyclic()} is a the column-wise analogue of
#' \code{as.rowcyclic()}.  In other words, it block-cyclically redistributes
#' the data across context 1.
#' 
#' \code{as.blockcyclic()} moves the distributed matrix into a general
#' block-cyclic distribution across a 2-dimensional process grid.  In other
#' words, it block-cyclically redistributes the data across context 0.
#' 
#' @param dx 
#' numeric distributed matrix
#' @param square.bldim 
#' logical.  Determines whether or not the blocking factor
#' for the resulting redistributed matrix will be square or not.
#' @param bldim 
#' the blocking dimension for block-cyclically distributing the
#' matrix across the process grid.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @keywords BLACS Distributing Data
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' dx <- ddmatrix(1:30, nrow=10)
#' 
#' x <- as.block(dx)
#' x
#' 
#' x <- as.rowblock(dx)
#' x
#' 
#' x <- as.colblock(dx)
#' x
#' 
#' x <- as.rowcyclic(dx)
#' x
#' 
#' x <- as.colcyclic(dx)
#' x
#' 
#' x <- as.blockcyclic(dx)
#' x
#' 
#' finalize()
#' }
#' 
#' @name as.rowcyclic
#' @rdname as.rowcyclic
#' @export
as.rowcyclic <- function(dx, bldim=.pbd_env$BLDIM)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  if (dx@ICTXT == 2 && all(dx@bldim == bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=2L)
  
  return( ret )
}



#' @rdname as.rowcyclic
#' @export
as.colcyclic <- function(dx, bldim=.pbd_env$BLDIM)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  if (dx@ICTXT == 1 && all(dx@bldim == bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=1L)
  
  return( ret )
}



#' @rdname as.rowcyclic
#' @export
as.blockcyclic <- function(dx, bldim=.pbd_env$BLDIM)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  if (dx@ICTXT == 0L && all(dx@bldim == bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=bldim, ICTXT=0L)
}



#' @rdname as.rowcyclic
#' @export
as.block <- function(dx, square.bldim=TRUE)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  blacs_ <- base.blacs(ICTXT=dx@ICTXT)
  procs <- c(blacs_$NPROW, blacs_$NPCOL)
  
  if (square.bldim)
    new.bldim <- rep(max(sapply(1L:2L, function(i) ceiling(dx@dim[i]/procs[i]))), 2L)
  else
    new.bldim <- sapply(1L:2L, function(i) ceiling(dx@dim[i]/procs[i]))
  
  if (all(dx@bldim == new.bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=new.bldim, ICTXT=dx@ICTXT)
}



#' @rdname as.rowcyclic
#' @export
as.rowblock <- function(dx)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  new.bldim <- rep(ceiling(dx@dim[1L]/comm.size()), 2L)
  
  if (dx@ICTXT == 2 && all(dx@bldim == new.bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=new.bldim, ICTXT=2L)
}



#' @rdname as.rowcyclic
#' @export
as.colblock <- function(dx)
{
  if (!is.ddmatrix(dx))
    stop("function only applies to objects of class 'ddmatrix'")
  
  new.bldim <- rep(ceiling(dx@dim[2L]/comm.size()), 2L)
  
  if (dx@ICTXT == 1 && all(dx@bldim == new.bldim))
    return(dx)
  else
    ret <- dmat.reblock(dx=dx, bldim=new.bldim, ICTXT=1L)
}
