#' GBD Matrix to Distributed Dense Matrix and vice versa
#' 
#' This function convert a GBD matrix and a distributed dense matrix.
#' 
#' \code{X.gbd} is a matrix with dimension \code{N.gbd * p} and exists on all
#' processors. \code{N.gbd} may be vary across processors.
#' 
#' If \code{skip.balance = TRUE}, then \code{load.balance} will not be called
#' and \code{X.gbd} is preassumed to be balanced.
#' 
#' For demonstration purpose, these objects should not contains weird values
#' such as \code{NA}.
#' 
#' \code{dmat2gbd} is supposed returned a balanced gbd matrix if
#' \code{bal.info} is not supplied.
#' 
#' @param X.gbd 
#' a GBD matrix.
#' @param skip.balance 
#' if \code{load.balance} were skipped.
#' @param comm 
#' a communicator number.
#' @param bldim 
#' the blocking dimension for block-cyclically distributing the
#' matrix across the process grid.
#' @param gbd.major 
#' 1 for row-major storage, 2 for column-major.
#' @param ICTXT 
#' BLACS context number for return.
#' @param X.dmat 
#' a ddmatrix matrix.
#' @param bal.info 
#' a returned object from \code{balance.info}.
#' 
#' @return 
#' \code{gbd2dmat} returns a ddmatrix object.  \code{dmat2gbd} returns
#' a (balanced) gbd matrix.
#' 
#' @examples
#' \dontrun{
#' ### Under command mode, run the demo with 4 processors by
#' ### (Use Rscript.exe for windows system)
#' mpiexec -np 4 Rscript -e "demo(gbd_dmat,'pbdDEMO',ask=F,echo=F)"
#' }
#' 
#' @keywords programming
#' @name gbd_dmat
#' @rdname gbd_dmat
NULL

demo.gbdr2dmat <- function(X.gbd, skip.balance = FALSE, comm = .pbd_env$SPMD.CT$comm,
    bldim = .pbd_env$BLDIM, ICTXT = .pbd_env$ICTXT){
  ### check data.
  if(! is.matrix(X.gbd)){
    dim(X.gbd) <- c(length(X.gbd), 1)
  }
  p <- as.integer(ncol(X.gbd))
  all.check <- spmd.allgather.integer(p, integer(comm.size(comm)), comm = comm)
  if(any(all.check != p)){
    comm.print(all.check)
    stop("X.gbd is not consistent accross processors.")
  }

  ### load balance data.
  if(! skip.balance){
    X.gbd <- load.balance(X.gbd, comm = comm, gbd.major = 1) 
  }
  ldim.gbd <- as.integer(dim(X.gbd))
  N <- spmd.allreduce.integer(ldim.gbd[1], integer(1), op = "sum", comm = comm)
  bldim.org <- c(spmd.allreduce.integer(ldim.gbd[1], integer(1), op = "max",
                                        comm = comm),
                 p)

  ### block-cyclic in context 2.
  ldim <- base.numroc(dim = c(N, p), bldim = bldim.org, ICTXT = 2)
  X.dmat <- new("ddmatrix", Data = X.gbd,
                dim = c(N, p), ldim = ldim, bldim = bldim.org, ICTXT = 2)

  ### reblock to any context and block size.
  X.dmat <- pbdDMAT::reblock(X.dmat, bldim = bldim, ICTXT = ICTXT)

  X.dmat
} # End of demo.gbdr2dmat().

demo.gbdc2dmat <- function(X.gbd, skip.balance = FALSE, comm = .pbd_env$SPMD.CT$comm,
    bldim = .pbd_env$BLDIM, ICTXT = .pbd_env$ICTXT){
  ### check data.
  if(! is.matrix(X.gbd)){
    dim(X.gbd) <- c(1, length(X.gbd))
  }
  p <- as.integer(nrow(X.gbd))
  all.check <- spmd.allgather.integer(p, integer(comm.size(comm)), comm = comm)
  if(any(all.check != p)){
    comm.print(all.check)
    stop("X.gbd is not consistent accross processors.")
  }

  ### load balance data.
  if(! skip.balance){
    X.gbd <- load.balance(X.gbd, comm = comm, gbd.major = 2) 
  }
  ldim.gbd <- as.integer(dim(X.gbd))
  N <- spmd.allreduce.integer(ldim.gbd[2], integer(1), op = "sum", comm = comm)
  bldim.org <- c(p, spmd.allreduce.integer(ldim.gbd[2], integer(1), op = "max",
                                           comm = comm))

  ### block-cyclic in context 1.
  ldim <- base.numroc(dim = c(p, N), bldim = bldim.org, ICTXT = 1)
  X.dmat <- new("ddmatrix", Data = X.gbd,
                dim = c(p, N), ldim = ldim, bldim = bldim.org, ICTXT = 1)

  ### reblock to any context and block size.
  X.dmat <- pbdDMAT::reblock(X.dmat, bldim = bldim, ICTXT = ICTXT)

  X.dmat
} # End of demo.gbdc2dmat().


#' @rdname gbd_dmat
#' @export
gbd2dmat <- function(X.gbd, skip.balance = FALSE, comm = .pbd_env$SPMD.CT$comm,
    gbd.major = .pbd_env$gbd.major, bldim = .pbd_env$BLDIM,
    ICTXT = .pbd_env$ICTXT){
  if(gbd.major == 1){
    demo.gbdr2dmat(X.gbd, skip.balance = skip.balance, comm = comm,
                    bldim = bldim, ICTXT = ICTXT)
  } else if(gbd.major == 2){
    demo.gbdc2dmat(X.gbd, skip.balance = skip.balance, comm = comm,
                    bldim = bldim, ICTXT = ICTXT)
  } else{
    stop("gbd.major = 1 or 2.")
  }
} # End of gbd2dmat().


demo.dmat2gbdr <- function(X.dmat, bal.info = NULL, comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- comm.size(comm)

  ### check data.
  all.check <- spmd.allreduce.integer(is.ddmatrix(X.dmat), integer(1),
                   op = "sum", comm = comm) == COMM.SIZE 
  if(!all.check){
    stop("X.dmat is not consistent accross processors.")
  }

  ### block-cyclic in context 2.
  bldim.new <- c(ceiling(nrow(X.dmat) / COMM.SIZE), ncol(X.dmat))
  X.dmat <- pbdDMAT::reblock(X.dmat, bldim = bldim.new, ICTXT = 2)

  ### copy to gbd.
  if(base.ownany(dim(X.dmat), bldim(X.dmat), ICTXT = 2)){
    X.gbd <- X.dmat@Data
  } else{
    X.gbd <- matrix(0, nrow = 0, ncol = 0)
  }

  ### unload balance data.
  if(! is.null(bal.info)){
    if(bal.info$gbd.major != 1){
      stop("gbd.major should be 1.")
    } 
    X.gbd <- unload.balance(X.gbd, bal.info, comm = comm)
  }

  X.gbd
} # End of demo.dmat2gbdr().

demo.dmat2gbdc <- function(X.dmat, bal.info = NULL, comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- comm.size(comm)

  ### check data.
  all.check <- spmd.allreduce.integer(is.ddmatrix(X.dmat), integer(1),
                   op = "sum", comm = comm) == COMM.SIZE 
  if(!all.check){
    stop("X.dmat is not consistent accross processors.")
  }

  ### block-cyclic in context 1.
  bldim.new <- c(nrow(X.dmat), ceiling(ncol(X.dmat) / COMM.SIZE))
  X.dmat <- pbdDMAT::reblock(X.dmat, bldim = bldim.new, ICTXT = 1)

  ### copy to gbd.
  if(base.ownany(dim(X.dmat), bldim(X.dmat), ICTXT = 1)){
    X.gbd <- X.dmat@Data
  } else{
    X.gbd <- matrix(0, nrow = 0, ncol = 0)
  }

  ### unload balance data.
  if(! is.null(bal.info)){
    if(bal.info$gbd.major != 2){
      stop("gbd.major should be 2.")
    } 
    X.gbd <- unload.balance(X.gbd, bal.info, comm = comm)
  }

  X.gbd
} # End of demo.dmat2gbdc().


#' @rdname gbd_dmat
#' @export
dmat2gbd <- function(X.dmat, bal.info = NULL, comm = .pbd_env$SPMD.CT$comm,
    gbd.major = .pbd_env$gbd.major){
  if(gbd.major == 1){
    demo.dmat2gbdr(X.dmat, bal.info = bal.info, comm = comm)
  } else if(gbd.major == 2){
    demo.dmat2gbdc(X.dmat, bal.info = bal.info, comm = comm)
  } else{
    stop("gbd.major = 1 or 2.")
  }
} # End of dmat2gbd().
