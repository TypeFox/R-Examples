### This is an independent function to provide logL.
### Note that Inf, -Inf, NA, NaN is drop from the summation.

indep.logL <- function(PARAM){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  nrow <- nrow(X.spmd)
  ncol <- ncol(X.spmd)

  log.a <- log(2 * pi) * (-0.5 * ncol)

  ret <- matrix(0, nrow = nrow, ncol = PARAM$K)
  for(i.k in 1:PARAM$K){
    tmp.X <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow, ncol)

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

  ret <- sum(log(ret[is.finite(ret)]))
  ret <- spmd.allreduce.double(ret, double(1), op = "sum")
  ret
} # End of indep.logL().

