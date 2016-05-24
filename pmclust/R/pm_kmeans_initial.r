### This file gives a simple initialization.

initial.center.spmd <- function(PARAM, MU = NULL){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  if(is.null(MU)){
    N.spmd <- nrow(X.spmd)
    N.allspmds <- spmd.allgather.integer(as.integer(N.spmd),
                                         integer(.pmclustEnv$COMM.SIZE))

    center.spmd <- rep(0, PARAM$K)
    if(.pmclustEnv$COMM.RANK == 0){
      center.spmd <- sample(1:.pmclustEnv$COMM.SIZE, PARAM$K, replace = TRUE,
                              prob = N.allspmds / PARAM$N) - 1
    }
    center.spmd <- spmd.bcast.integer(as.integer(center.spmd))

    tmp <- NULL
    n.center.spmd <- sum(center.spmd == .pmclustEnv$COMM.RANK)
    if(n.center.spmd > 0){
      id.center.spmd <- sample(1:N.spmd, n.center.spmd)
      tmp <- matrix(X.spmd[id.center.spmd,], ncol = ncol(X.spmd),
                    byrow = TRUE)
    }

    PARAM$MU <- unlist(spmd.allgather.object(tmp))
    dim(PARAM$MU) <- c(PARAM$p, PARAM$K)
  } else{
    PARAM$MU <- MU
  }

  for(i.k in 1:PARAM$K){
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow(X.spmd), ncol(X.spmd))
    .pmclustEnv$Z.spmd[, i.k] <- -rowSums(B * B)
  }

  .pmclustEnv$CLASS.spmd <- unlist(apply(.pmclustEnv$Z.spmd, 1, which.max))

  PARAM
} # End of initial.center.spmd().

