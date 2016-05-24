### For general internal methods.

pmclust.internal.dmat <- function(X = NULL, K = 2, MU = NULL,
    algorithm = .PMC.CT$algorithm.dmat, RndEM.iter = .PMC.CT$RndEM.iter,
    CONTROL = .PMC.CT$CONTROL, method.own.X = .PMC.CT$method.own.X,
    rank.own.X = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  # Check.
  if(! (algorithm[1] %in% .PMC.CT$algorithm.dmat)){
    comm.stop("The algorithm is not supported.")
  }

  # Check X.
  if(comm.all(is.null(X))){
    if(! pbdDMAT::is.ddmatrix(.GlobalEnv$X.dmat)){
      comm.stop("X.dmat is not a ddmatrix.")
    } else{
      # Assume X.dmat in .GlobalEnv and no need for converting.
      PARAM.org <- set.global.dmat(K = K, RndEM.iter = RndEM.iter)
    }
  } else{
    # Assign X to .pmclustEnv if it is not in .GlobalEnv
    if(! pbdDMAT::is.ddmatrix(X)){
      comm.stop("A ddmatrix is required.")
    }
    PARAM.org <- set.global.dmat(K = K, RndEM.iter = RndEM.iter, X.dmat = X)
  }

  # Set global variables.
  # PARAM.org <- set.global.dmat(K = K, RndEM.iter = RndEM.iter)
  if(! comm.all(is.null(CONTROL))){
    tmp <- .pmclustEnv$CONTROL[!(names(.pmclustEnv$CONTROL) %in%
                                 names(CONTROL))]
    .pmclustEnv$CONTROL <- c(tmp, CONTROL)
  }

  # Initialization for algorithms.
  if(! comm.all(is.null(MU))){
    if(algorithm[1] != "kmeans.dmat"){
      PARAM.org <- initial.em.dmat(PARAM.org, MU = MU)
    } else{
      PARAM.org <- initial.center.dmat(PARAM.org, MU = MU)
    }
  } else{
    if(algorithm[1] != "kmeans.dmat"){
      PARAM.org <- initial.RndEM.dmat(PARAM.org)
    } else{
      PARAM.org <- initial.center.dmat(PARAM.org)
    }
  }

  # Update steps.
  method.step <- switch(algorithm[1],
                        "em.dmat" = em.step.dmat,
#                        "aecm.dmat" = aecm.step.dmat,
#                        "apecm.dmat" = apecm.step.dmat,
#                        "apecma.dmat" = apecma.step.dmat,
                        "kmeans.dmat" = kmeans.step.dmat,
                        NULL)
  if(comm.all(is.null(method.step))){
    comm.stop("Algorithm is not found.")
  }
  PARAM.new <- method.step(PARAM.org)

  # Obtain classifications.
  if(algorithm[1] == "kmeans"){
    kmeans.update.class.dmat()
  } else{
    em.update.class.dmat()
  }

  # Get class numbers.
  N.CLASS <- get.N.CLASS.dmat(K)

  # For return.
  ret <- list(algorithm = algorithm[1],
              param = PARAM.new,
              class = .pmclustEnv$CLASS.spmd,
              n.class = N.CLASS,
              check = .pmclustEnv$CHECK)
  ret
} # End of pmclust.internal().

