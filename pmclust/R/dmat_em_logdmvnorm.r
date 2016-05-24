### This file contains functions for log density of MVN.
### These will majorly update .pmclustEnv$W.spmd.

logdmvnorm.dmat <- function(PARAM, i.k){
  if(exists("X.dmat", envir = .pmclustEnv)){
    X.dmat <- get("X.dmat", envir = .pmclustEnv)
  }

#  for(i.k in 1:PARAM$K){
#    U <- chol(PARAM$SIGMA[[i.k]])
    U <- PARAM$U[[i.k]]
    logdet <- sum(log(abs(diag(U)))) * 2
#    B <- t.X.spmd - PARAM$MU[, i.k]
#    A <- backsolve(U, B, upper.tri = TRUE, transpose = TRUE)
#    distval <- colSums(A * A)

    ### WCC: original
    # B <- sweep(X.dmat, 2, PARAM$MU[, i.k])
    # C <- backsolve(U, diag(1, PARAM$p))
    # B <- B %*% pbdDMAT::as.ddmatrix(C, bldim = bldim(X.dmat), ICTXT = pbdDMAT::ICTXT(X.dmat))
    # distval <- rowSums(B * B)
    # .pmclustEnv$W.dmat[, i.k] <- -(.pmclustEnv$p.times.logtwopi + logdet +
    #                                distval) * 0.5
    ### WCC: temp dmat
    # tmp.1 <- sweep(X.dmat, 2, PARAM$MU[, i.k])
    # tmp.2 <- backsolve(U, diag(1, PARAM$p))
    # tmp.3 <- pbdDMAT::as.ddmatrix(tmp.2, bldim = bldim(X.dmat), ICTXT = pbdDMAT::ICTXT(X.dmat))
    # tmp.4 <- tmp.1 %*% tmp.3
    # tmp.5 <- tmp.4 * tmp.4
    # tmp.6 <- rowSums(tmp.5)
    # tmp.7 <- -(.pmclustEnv$p.times.logtwopi + logdet + tmp.6) * 0.5
    # .pmclustEnv$W.dmat[, i.k] <- tmp.7
    ### WCC: temp spmd
    tmp.1 <- as.matrix(X.dmat)
    B <- sweep(tmp.1, 2, PARAM$MU[, i.k])
    B <- B %*% backsolve(U, diag(1, PARAM$p))
    distval <- rowSums(B * B)
    tmp.2 <- -(.pmclustEnv$p.times.logtwopi + logdet + distval) * 0.5
    tmp.3 <- as.matrix(.pmclustEnv$W.dmat)
    tmp.3[, i.k] <- tmp.2
    .pmclustEnv$W.dmat <- pbdDMAT::as.ddmatrix(tmp.3)

    # B <- sweep(X.dmat, 2, PARAM$MU[, i.k])
    # C <- backsolve(U, diag(1, PARAM$p))
    # B <- B %*% pbdDMAT::as.ddmatrix(C)
    # distval <- rowSums(B * B)
    # distval <- as.vector(distval)
    # .pmclustEnv$W.dmat[, i.k] <- -(.pmclustEnv$p.times.logtwopi + logdet +
    #                                distval) * 0.5

#  }
  invisible()
} # End of logdmvnorm.dmat().

