### This function initializes global variables.
set.global.dmat <- function(K = 2, X.dmat = NULL, PARAM = NULL,
    algorithm = c("em.dmat", "kmeans.dmat"),
    RndEM.iter = 10){
  if(is.null(X.dmat)){
    if(exists("X.dmat", envir = .pmclustEnv)){
      X.dmat <- get("X.dmat", envir = .pmclustEnv)
    } else{
      if(! exists("X.dmat", envir = .GlobalEnv)){
        comm.stop("A global X.dmat does not exist.")
      }
    }
  } else{
    .pmclustEnv$X.dmat <- X.dmat
  }

  if(! pbdDMAT::is.ddmatrix(X.dmat)){
    stop("X.dmat is not a ddmatrix.")
  }
  CTXT <- pbdDMAT::ICTXT(X.dmat)

  ### Get data information.
  N <- nrow(X.dmat)
  p <- ncol(X.dmat)

  ### Set parameters.
  if(is.null(PARAM)){
    PARAM <- list(N = N, p = p, K = K,
                  ETA = NULL, log.ETA = NULL, MU = NULL, SIGMA = NULL,
                  U = rep(list(), K),
                  U.check = rep(TRUE, K),
                  logL = NULL,
                  min.N.CLASS = min(c((p + 1) * p * 0.5 + 1, N / K * 0.2)))
    PARAM$ETA <- rep(1 / K, K)
    PARAM$log.ETA <- rep(-log(K), K) 
    PARAM$MU <- matrix(0, p, K)
    PARAM$SIGMA <- rep(list(diag(1.0, p)), K)
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

  .pmclustEnv$Z.dmat <- pbdDMAT::ddmatrix(0, N, K)
  .pmclustEnv$Z.colSums <- as.vector(colSums(.pmclustEnv$Z.dmat))

  .pmclustEnv$W.dmat <- pbdDMAT::ddmatrix(0, N, K)
  .pmclustEnv$W.rowSums <- as.vector(rowSums(.pmclustEnv$W.dmat))

  .pmclustEnv$U.dmat <- pbdDMAT::ddmatrix(0, N, K)

  # .pmclustEnv$CLASS.dmat <- pbdDMAT::ddmatrix(0, N, 1)
  .pmclustEnv$CLASS <- rep(0, N)	# This is not a ddmatrix.

  .pmclustEnv$CHECK <- list(algorithm = algorithm[1],
                            i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)

  ### For semi-supervised clustering.
#  assign.ss.spmd()

  for(i.k in 1:K){
    tmp.U <- decompsigma(PARAM$SIGMA[[i.k]])
    PARAM$U[[i.k]] <- tmp.U$value 
    PARAM$U.check[[i.k]] <- tmp.U$check 
  }

  PARAM
} # End of set.global.dmat().

