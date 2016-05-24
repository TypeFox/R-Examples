### This is an independent function to provide logL.
### Note that Inf, -Inf, NA, NaN is drop from the summation.

indep.logL.dmat <- function(PARAM){
  if(exists("X.dmat", envir = .pmclustEnv)){
    X.dmat <- get("X.dmat", envir = .pmclustEnv)
  }

  nrow <- nrow(X.dmat)
  ncol <- ncol(X.dmat)

  log.a <- log(2 * pi) * (-0.5 * ncol)

  ret <- matrix(0, nrow = nrow, ncol = PARAM$K)
  for(i.k in 1:PARAM$K){
    ### bug
    #tmp.X <- base.pdsweep(dx = X.dmat, vec = PARAM$MU[, i.k],
    #                      MARGIN = 2L, FUN = "-")
    comm.stop("Not implemented yet.")
    B <- NULL
    tmp.X <- as.matrix(tmp.X)

    tmp.S <- PARAM$SIGMA[[i.k]]
    log.b <- -0.5 * log(abs(det(tmp.S)))

    tmp.S <- solve(tmp.S)
    log.c <- -0.5 * rowSums((tmp.X %*% tmp.S) * tmp.X)

    ret[, i.k] <- log.c + log.b + log.a + PARAM$log.ETA[i.k]
  }

  ret <- rowSums(exp(ret))

  if(.pmclustEnv$CONTROL$debug > 10){
    comm.cat("  >>Not finite: ", sep = "", quiet = TRUE)
    for(i.rank in 0:(.pmclustEnv$COMM.SIZE - 1)){
      if(i.rank == .pmclustEnv$COMM.RANK){
        cat(.pmclustEnv$COMM.RANK, ":", sum(!is.finite(ret)), " ", sep = "")
      }
      # spmd.barrier()
    }
    comm.cat("\n", sep = "", quiet = TRUE)
  }

  sum(log(ret[is.finite(ret)]))
} # End of indep.logL.dmat().

