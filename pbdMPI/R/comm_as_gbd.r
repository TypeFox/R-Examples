### Redistributing a regular matrix to be a gbd matrix.
### Assume gbd.major = 1.

comm.as.gbd <- function(X, balance.method = .pbd_env$SPMD.IO$balance.method,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  ### Check.
  COMM.RANK <- spmd.comm.rank(comm)
  ncol.X <- integer(1)
  if(COMM.RANK == rank.source){
    if(!is.matrix(X)){
      stop("X should be a matrix.")
    }
    ncol.X <- ncol(X)
  }
  ncol.X <- spmd.bcast.integer(ncol.X, rank.source = rank.source, comm = comm)

  ### Elimate all others.
  if(COMM.RANK != rank.source){
    X <- matrix(0, nrow = 0, ncol = ncol.X)
  }

  ### Redistributed.
  ret <- comm.load.balance(X, balance.method = balance.method, comm = comm)
  ret
} # End of comm.as.gdb().
