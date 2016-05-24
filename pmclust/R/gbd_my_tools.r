### This function initializes global variables.
set.global.gbd <- function(K = 2, X.gbd = NULL, PARAM = NULL,
    algorithm = c("em", "aecm", "apecm", "apecma", "kmeans"),
    RndEM.iter = 10){
  if(is.null(X.gbd)){
    if(exists("X.gbd", envir = .pmclustEnv)){
      .pmclustEnv$X.spmd <- X.gbd
    } else{
      A <- exists("X.gbd", envir = .GlobalEnv)
      B <- exists("X.spmd", envir = .GlobalEnv)
      if(A){
        .pmclustEnv$X.spmd <- get("X.gbd", envir = .GlobalEnv)
      } else if(B){
        .pmclustEnv$X.spmd <- get("X.spmd", envir = .GlobalEnv)
      } else{
        comm.stop("A global X.gbd or X.spmd does not exist.")
      }
    }
  } else{
    .pmclustEnv$X.spmd <- X.gbd
  }
  X.spmd <- .pmclustEnv$X.spmd

  ### Get data information.
  N.spmd <- nrow(X.spmd)
  N.allspmds <- spmd.allgather.integer(as.integer(N.spmd),
                                       integer(spmd.comm.size()))
  N <- sum(N.allspmds)
  p <- ncol(X.spmd)

  ### Set parameters.
  if(is.null(PARAM)){
    PARAM <- list(N = N, p = p, K = K, ETA = rep(1 / K, K),
                  log.ETA = rep(-log(K), K), MU = NULL,
                  SIGMA = rep(list(diag(1.0, p)), K),
                  U = rep(list(), K),
                  U.check = rep(TRUE, K),
                  logL = NULL,
                  min.N.CLASS = min(c((p + 1) * p * 0.5 + 1, N / K * 0.2)))
  } else{
    PARAM$N <- N
    K <- PARAM$K
  }

  ### Set global storages.
  .pmclustEnv$CONTROL <- .PMC.CT$CONTROL
  .pmclustEnv$CONTROL$RndEM.iter <- RndEM.iter

  .pmclustEnv$COMM.SIZE <- spmd.comm.size()
  .pmclustEnv$COMM.RANK <- spmd.comm.rank()

  .pmclustEnv$p.times.logtwopi <- p * log(2 * pi)

  .pmclustEnv$Z.spmd <- matrix(0.0, nrow = N.spmd, ncol = K)
  .pmclustEnv$Z.colSums <- rep(0.0, K)

  .pmclustEnv$W.spmd <- matrix(0.0, nrow = N.spmd, ncol = K)
  .pmclustEnv$W.spmd.rowSums <- rep(0.0, N.spmd)

  .pmclustEnv$U.spmd <- matrix(0.0, nrow = N.spmd, ncol = K)
  .pmclustEnv$CLASS.spmd <- rep(0, N.spmd)
  .pmclustEnv$CHECK <- list(algorithm = algorithm[1], i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)

  ### For semi-supervised clustering.
#  assign.ss.spmd()

  for(i.k in 1:K){
    tmp.U <- decompsigma(PARAM$SIGMA[[i.k]])
    PARAM$U[[i.k]] <- tmp.U$value 
    PARAM$U.check[[i.k]] <- tmp.U$check 
  }

  PARAM
} # End of set.global.gbd().
