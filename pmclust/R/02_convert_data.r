### For general methods.

convert.data <- function(X, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- spmd.comm.size(comm)
  COMM.RANK <- spmd.comm.rank(comm)

  # Assign X to .pmclustEnv
  if(pbdDMAT::is.ddmatrix(X)){
    # For a ddmatrix.

    .pmclustEnv$X.spmd <- as.spmd(X, comm = comm)
  } else{
    # for a spmd matrix.

    if(is.null(X)){
      X <- matrix(0, nrow = 0, ncol = 0)
    }

    if(method.own.X[1] %in% c("gbdr", "spmdr")){
      # For spmd row-major

      p <- ncol(X)
      p.all <- spmd.allgather.integer(p, integer(COMM.SIZE), comm = comm)
      if(any(p.all != p)){
        comm.stop("X should have the same # of columns.")
      }

      .pmclustEnv$X.spmd <- X
    } else if(method.own.X[1] == "common"){
      # X is common in all ranks.

      N <- nrow(X)
      N.all <- spmd.allgather.integer(N, integer(COMM.SIZE), comm = comm)
      if(any(N.all != N)){
        comm.stop("X should have the same # of rows.")
      }

      p <- ncol(X)
      p.all <- spmd.allgather.integer(p, integer(COMM.SIZE), comm = comm)
      if(any(p.all != p)){
        comm.stop("X should have the same # of columns.")
      }

      jid <- get.jid(nrow(X))
      .pmclustEnv$X.spmd <- X[jid,]
    } else if(method.own.X[1] == "single"){
      # X only exist in a single rank, 0 by default.

      p <- -1
      if(COMM.RANK == rank.own.X){
        if(is.matrix(X)){
          p <- ncol(X)
        }
      }

      p <- spmd.bcast.integer(as.integer(p), rank.source = rank.own.X,
                              comm = comm)
      if(p == -1){
        comm.stop("X should be a matrix in rank.own.X.")
      } else{
        if(COMM.RANK != rank.own.X){
          X <- matrix(0, nrow = 0, ncol = p)
        }
      }

      .pmclustEnv$X.spmd <- load.balance(X)
    } else{
      comm.stop("method.own.X is not found.")
    }
  }

  invisible()
} # End of convert.data().
