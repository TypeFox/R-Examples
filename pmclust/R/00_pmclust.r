### For general methods.

pmclust <- function(X = NULL, K = 2, MU = NULL,
    algorithm = .PMC.CT$algorithm, RndEM.iter = .PMC.CT$RndEM.iter,
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  if(comm.all(is.null(X))){
    # Check global matrix.
    A <- exists("X.spmd", envir = .GlobalEnv)
    B <- exists("X.dmat", envir = .GlobalEnv)
    if((!A) & (!B)){
      if(! algorithm[1] %in% .PMC.CT$algorithm.gbd){
        comm.stop("A global X.spmd is required in .GlobalEnv.")
      }
      if(! algorithm[1] %in% .PMC.CT$algorithm.dmat){
        comm.stop("A global X.dmat is required in .GlobalEnv.")
      }
    }
    if(A & B){
      comm.stop("Both X.spmd and X.dmat are in .GlobalEnv.")
    }
    # Check matrix type if dmat algorithm is used.
    if(B & algorithm[1] %in% .PMC.CT$algorithm.dmat){
      if(! pbdDMAT::is.ddmatrix(.GlobalEnv$X.dmat)){
        comm.stop("X.dmat is not a ddmatrix.")
      }
    }
  } else{
    # Check matrix type if dmat algorithm is used.
    if(algorithm[1] %in% .PMC.CT$algorithm.dmat){
      if(! pbdDMAT::is.ddmatrix(X)){
        comm.stop("X is not a ddmatrix.")
      }
    }
  }

  if(algorithm[1] %in% .PMC.CT$algorithm.gbd){
    ret <- pmclust.internal(X, K,
                            MU = MU,
                            algorithm = algorithm[1],
                            RndEM.iter = RndEM.iter,
                            CONTROL = CONTROL,
                            method.own.X = method.own.X[1],
                            rank.own.X = rank.own.X,
                            comm = comm)
  } else if(algorithm[1] %in% .PMC.CT$algorithm.dmat){
    ret <- pmclust.internal.dmat(X, K,
                                 MU = MU,
                                 algorithm = algorithm[1],
                                 RndEM.iter = RndEM.iter,
                                 CONTROL = CONTROL,
                                 method.own.X = method.own.X[1],
                                 rank.own.X = rank.own.X,
                                 comm = comm)
  } else{
    comm.stop("The algorithm is not found.")
  }

  if(algorithm[1] %in% c("kmeans", "kmeans.dmat")){
    class(ret) <- "pkmeans"
  } else{
    class(ret) <- "pmclust"
  }
  ret
} # end of pmclust().


pkmeans <- function(X = NULL, K = 2, MU = NULL,
    algorithm = c("kmeans", "kmeans.dmat"),
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  if(comm.all(is.null(X))){
    # Check global matrix.
    A <- exists("X.spmd", envir = .GlobalEnv)
    B <- exists("X.dmat", envir = .GlobalEnv)
    if((!A) & (!B)){
      if(! algorithm[1] %in% .PMC.CT$algorithm.gbd){
        comm.stop("A global X.spmd is required in .GlobalEnv.")
      }
      if(! algorithm[1] %in% .PMC.CT$algorithm.dmat){
        comm.stop("A global X.dmat is required in .GlobalEnv.")
      }
    }
    if(A & B){
      comm.stop("Both X.spmd and X.dmat are in .GlobalEnv.")
    }
    # Check matrix type if dmat algorithm is used.
    if(B & algorithm[1] %in% .PMC.CT$algorithm.dmat){
      if(! pbdDMAT::is.ddmatrix(.GlobalEnv$X.dmat)){
        comm.stop("X.dmat is not a ddmatrix.")
      }
    }
  } else{
    # Check matrix type if dmat algorithm is used.
    if(algorithm[1] %in% .PMC.CT$algorithm.dmat){
      if(! pbdDMAT::is.ddmatrix(X)){
        comm.stop("X is not a ddmatrix.")
      }
    }
  }

  if(algorithm[1] == "kmeans"){
    ret <- pmclust.internal(X, K,
                            MU = MU,
                            algorithm = algorithm[1],
                            CONTROL = CONTROL,
                            method.own.X = method.own.X[1],
                            rank.own.X = rank.own.X,
                            comm = comm)
  } else if(algorithm[1] == "kmeans.dmat"){
    ret <- pmclust.internal.dmat(X, K,
                                 MU = MU,
                                 algorithm = algorithm[1],
                                 CONTROL = CONTROL,
                                 method.own.X = method.own.X[1],
                                 rank.own.X = rank.own.X,
                                 comm = comm)
  } else{
    comm.stop("The algorithm is not found.")
  }
  class(ret) <- "pkmeans"
  ret
} # end of pkmeans().

