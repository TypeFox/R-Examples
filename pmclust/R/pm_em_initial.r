### This file gives initializations.

initial.em.spmd <- function(PARAM, MU = NULL){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  if(is.null(MU)){
    ### Set semi-supervised information.
    N.spmd <- nrow(X.spmd)
    unlabeled.K <- PARAM$K
    unlabeled.N.spmd <- N.spmd
    # if(SS.clustering){
    #   unlabeled.K <- unlabeled.K - SS.K
    #   unlabeled.N.spmd <- unlabeled.N.spmd - length(SS.id.spmd)
    # }

    ### Simple random sampling from data.
    N.allspmds <- spmd.allgather.integer(as.integer(unlabeled.N.spmd),
                                         integer(.pmclustEnv$COMM.SIZE))

    center.spmd <- rep(0, unlabeled.K)
    if(.pmclustEnv$COMM.RANK == 0){
      center.spmd <- sample(1:.pmclustEnv$COMM.SIZE, unlabeled.K,
                            replace = TRUE,
                            prob = N.allspmds / sum(N.allspmds)) - 1
    }
    center.spmd <- spmd.bcast.integer(as.integer(center.spmd))

    tmp <- NULL
    n.center.spmd <- sum(center.spmd == .pmclustEnv$COMM.RANK)
    if(n.center.spmd > 0){
      N.pool <- 1:N.spmd
      # if(SS.clustering && length(SS.id.spmd) > 0){
      #   N.pool <- N.pool[-SS.id.spmd]
      # }
      id.center.spmd <- sample(N.pool, n.center.spmd)
      tmp <- matrix(X.spmd[id.center.spmd,], ncol = ncol(X.spmd),
                    byrow = TRUE)
    }
    MU <- matrix(unlist(spmd.allgather.object(tmp)),
                 nrow = ncol(X.spmd), ncol = unlabeled.K)

    ### Combind centers from semi-supervised information if any.
    # PARAM$MU <- cbind(SS.MU, MU)
    PARAM$MU <- MU
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
# PARAM$MU <- c(
#   0.2328901,  0.7281706, -1.0026075,
#  -0.6238820, -0.2584021,  0.8862744,
#   0.4944164,  0.7662080, -1.2993873,
#   0.4674663,  0.7461591, -1.2516524
# )
# PARAM$MU <- matrix(PARAM$MU, nrow = 4)

  e.step.spmd(PARAM)
  PARAM <- em.onestep.spmd(PARAM)
  PARAM$logL <- logL.step.spmd()
  em.update.class.spmd()

  PARAM
} # End of initial.em.spmd().

initial.RndEM.spmd <- function(PARAM){
  logL.save <- -Inf
  i.iter <- 1

  PARAM.org <- PARAM
  repeat{
    PARAM <- try(initial.em.spmd(PARAM.org))
    if(comm.any(class(PARAM) == "try-error")){
      comm.cat(PARAM, "\n", quiet = TRUE)
      next
    }

    N.CLASS <- get.N.CLASS(PARAM$K)
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
  PARAM <- initial.em.spmd(PARAM.save, MU = PARAM.save$MU)
  PARAM
} # End of initial.RndEM.spmd().

