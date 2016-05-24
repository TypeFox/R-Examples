### This file gives initializations.

initial.em.dmat <- function(PARAM, MU = NULL){
  if(exists("X.dmat", envir = .pmclustEnv)){
    X.dmat <- get("X.dmat", envir = .pmclustEnv)
  }

  if(! pbdDMAT::is.ddmatrix(X.dmat)){
    stop("X.dmat is not a ddmatrix.")
  }

  if(is.null(MU)){
    N <- nrow(X.dmat)
    id <- spmd.bcast.integer(as.integer(sample(1:N, PARAM$K)))
    PARAM$MU <- t(as.matrix(X.dmat[id, ]))
  } else{
    PARAM$MU <- MU
  }

### For iris example.
# PARAM$MU <- c(
# -0.8976739, 1.3968289, 0.5514857,
#  1.0156020, 0.3273175, 0.5567457,
# -1.3357516, 0.5336209, 1.2700404,
# -1.3110521, 0.2632600, 1.7063794
# )
#PARAM$MU <- c(
#  0.2328901,  0.7281706, -1.0026075,
# -0.6238820, -0.2584021,  0.8862744,
#  0.4944164,  0.7662080, -1.2993873,
#  0.4674663,  0.7461591, -1.2516524
#)
#PARAM$MU <- matrix(PARAM$MU, nrow = 4)

  e.step.dmat(PARAM)
  PARAM <- em.onestep.dmat(PARAM)
  PARAM$logL <- logL.step.dmat()
  em.update.class.dmat()

  PARAM
} # End of initial.em.dmat().

initial.RndEM.dmat <- function(PARAM){
  logL.save <- -Inf
  i.iter <- 1

  PARAM.org <- PARAM
  repeat{
    PARAM <- try(initial.em.dmat(PARAM.org))
    if(comm.any(class(PARAM) == "try-error")){
      comm.cat(PARAM, "\n", quiet = TRUE)
      next
    }

    N.CLASS <- get.N.CLASS.dmat(PARAM$K)
    if(any(N.CLASS < PARAM$min.N.CLASS)){
      if(.pmclustEnv$CONTROL$debug > 0){
        comm.cat("N.CLASS: ", N.CLASS, "\n", quiet = TRUE)
      }
      next
    }

    if(.pmclustEnv$CONTROL$debug > 0){
      comm.cat("Initial: ", format(Sys.time(), "%H:%M:%S"),
               ", iter: ", i.iter, ", logL: ",
                           sprintf("%-20.10f", PARAM$logL), "\n",
               sep = "", quiet = TRUE)
    }

    if(logL.save < PARAM$logL){
      logL.save <- PARAM$logL
      PARAM.save <- PARAM
      PARAM.save$initial.i.iter <- i.iter
    }

    i.iter <- i.iter + 1
    if(i.iter > .pmclustEnv$CONTROL$RndEM.iter){
      break
    }
  }

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat("Using initial iter: ", PARAM.save$initial.i.iter, "\n",
             sep = "", quiet = TRUE)
  }
  PARAM <- initial.em.dmat(PARAM.save, MU = PARAM.save$MU)
  PARAM
} # End of initial.RndEM.dmat().

