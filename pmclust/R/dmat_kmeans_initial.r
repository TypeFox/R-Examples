# This file gives a simple initialization.

initial.center.dmat <- function(PARAM, MU = NULL){
  if(exists("X.dmat", envir = .pmclustEnv)){
    X.dmat <- get("X.dmat", envir = .pmclustEnv)
  }

  if(! pbdDMAT::is.ddmatrix(X.dmat)){
    stop("X.dmat is not a ddmatrix.")
  }

  if(is.null(MU)){
    N <- nrow(X.dmat)
    id <- spmd.bcast.integer(as.integer(sample(1:N, PARAM$K)))
    ### WCC: original
    PARAM$MU <- t(as.matrix(X.dmat[id, ]))
    ### WCC: temp
    # tmp.1 <- X.dmat[id,]
    # tmp.2 <- as.matrix(tmp.1)
    # tmp.3 <- t(tmp.2)
    # PARAM$MU <- tmp.3 
  } else{
    PARAM$MU <- MU
  }

  for(i.k in 1:PARAM$K){
     ### WCC: original
     B <- sweep(X.dmat, 2, PARAM$MU[, i.k])
     .pmclustEnv$Z.dmat[, i.k] <- -rowSums(B * B)
     ### WCC: temp
     # tmp.1 <- sweep(X.dmat, 2, PARAM$MU[, i.k])
     # tmp.2 <- tmp.1 * tmp.1
     # tmp.3 <- rowSums(tmp.2)
     # tmp.4 <- -tmp.3
     # .pmclustEnv$Z.dmat[, i.k] <- tmp.4
  }

  ### WCC: original
  # .pmclustEnv$CLASS.dmat <- apply(.pmclustEnv$Z.dmat, 1, which.max)
  ### WCC: temp
  tmp.1 <- as.matrix(.pmclustEnv$Z.dmat)
  .pmclustEnv$CLASS <- unlist(apply(tmp.1, 1, which.max))

  PARAM
} # End of initial.center.dmat().

