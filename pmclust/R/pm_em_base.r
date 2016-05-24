### This file contains major functions for EM iterations.

### E-step.
e.step.spmd <- function(PARAM, update.logL = TRUE){
  for(i.k in 1:PARAM$K){
    logdmvnorm(PARAM, i.k)
  }

  update.expectation(PARAM, update.logL = update.logL)
  invisible()
} # End of e.step.spmd().

### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
update.expectation <- function(PARAM, update.logL = TRUE){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  N <- nrow(X.spmd)
  K <- PARAM$K

  .pmclustEnv$U.spmd <- W.plus.y(.pmclustEnv$W.spmd, PARAM$log.ETA, N, K)
  .pmclustEnv$Z.spmd <- exp(.pmclustEnv$U.spmd)

  tmp.id <- rowSums(.pmclustEnv$U.spmd < .pmclustEnv$CONTROL$exp.min) == K |
            rowSums(.pmclustEnv$U.spmd > .pmclustEnv$CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.spmd <- .pmclustEnv$U.spmd[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.spmd) - .pmclustEnv$CONTROL$exp.max / K
    } else{
      tmp.scale <- unlist(apply(tmp.spmd, 1, max)) -
                   .pmclustEnv$CONTROL$exp.max / K
    }
    .pmclustEnv$Z.spmd[tmp.id,] <- exp(tmp.spmd - tmp.scale)
  }

  .pmclustEnv$W.spmd.rowSums <- rowSums(.pmclustEnv$Z.spmd)
  .pmclustEnv$Z.spmd <- .pmclustEnv$Z.spmd / .pmclustEnv$W.spmd.rowSums

  ### For semi-supervised clustering.
  # if(SS.clustering){
  #   .pmclustEnv$Z.spmd[SS.id.spmd,] <- SS..pmclustEnv$Z.spmd
  # }

  .pmclustEnv$Z.colSums <- colSums(.pmclustEnv$Z.spmd)
  .pmclustEnv$Z.colSums <- spmd.allreduce.double(.pmclustEnv$Z.colSums,
                                                 double(K), op = "sum")

  if(update.logL){
    .pmclustEnv$W.spmd.rowSums <- log(.pmclustEnv$W.spmd.rowSums)
    if(tmp.flag > 0){
      .pmclustEnv$W.spmd.rowSums[tmp.id] <- .pmclustEnv$W.spmd.rowSums[tmp.id] +
                                            tmp.scale
    }
  }

  invisible()
} # End of update.expectation().


### M-step.
m.step.spmd <- function(PARAM){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  ### MLE For ETA
  PARAM$ETA <- .pmclustEnv$Z.colSums / sum(.pmclustEnv$Z.colSums)
  PARAM$log.ETA <- log(PARAM$ETA)

  p <- PARAM$p
  p.2 <- p * p
  for(i.k in 1:PARAM$K){
    ### MLE for MU
    tmp.MU <- colSums(X.spmd * .pmclustEnv$Z.spmd[, i.k]) /
              .pmclustEnv$Z.colSums[i.k]
    PARAM$MU[, i.k] <- spmd.allreduce.double(tmp.MU, double(p), op = "sum")

    ### MLE for SIGMA
    if(PARAM$U.check[[i.k]]){
      B <- W.plus.y(X.spmd, -PARAM$MU[, i.k],
                    nrow(X.spmd), ncol(X.spmd)) *
           sqrt(.pmclustEnv$Z.spmd[, i.k] / .pmclustEnv$Z.colSums[i.k])
      tmp.SIGMA <- crossprod(B)
      tmp.SIGMA <- spmd.allreduce.double(tmp.SIGMA, double(p.2), op = "sum") 
      dim(tmp.SIGMA) <- c(p, p)

      tmp.U <- decompsigma(tmp.SIGMA)
      PARAM$U.check[[i.k]] <- tmp.U$check
      if(tmp.U$check){
        PARAM$U[[i.k]] <- tmp.U$value
        PARAM$SIGMA[[i.k]] <- tmp.SIGMA
      }
    } else{
      if(.pmclustEnv$CONTROL$debug > 2){
        comm.cat("  SIGMA[[", i.k, "]] is fixed.\n", sep = "", quiet = TRUE)
      }
    }
  }

  PARAM
} # End of m.step.spmd().


### log likelihood.
logL.step.spmd <- function(){
  tmp.logL <- sum(.pmclustEnv$W.spmd.rowSums)
  spmd.allreduce.double(tmp.logL, double(1), op = "sum")
} # End of logL.step.spmd().


### Check log likelihood convergence.
check.em.convergence <- function(PARAM.org, PARAM.new, i.iter){
  abs.err <- PARAM.new$logL - PARAM.org$logL
  rel.err <- abs.err / abs(PARAM.org$logL)
  convergence <- 0

  if(abs.err < 0){
    convergence <- 4
  } else if(any(PARAM.new$ETA < PARAM.new$min.N.CLASS / PARAM.new$N)){
    convergence <- 3
  } else if(i.iter > .pmclustEnv$CONTROL$max.iter){
    convergence <- 2
  } else if(rel.err < .pmclustEnv$CONTROL$rel.err){
    convergence <- 1
  }

  if(.pmclustEnv$CONTROL$debug > 1){
    comm.cat("  check.em.convergence:",
             " abs: ", abs.err,
             ", rel: ", rel.err,
             ", conv: ", convergence, "\n",
             sep = "", quiet = TRUE)
  }

  list(algorithm = .pmclustEnv$CHECK$algorithm,
       iter = i.iter, abs.err = abs.err, rel.err = rel.err,
       convergence = convergence)
} # End of check.em.convergence().


### EM-step.
em.step.spmd <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(algorithm = "em", i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(.pmclustEnv$CONTROL$save.log)) && .pmclustEnv$CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .pmclustEnv)){
      .pmclustEnv$SAVE.param <- NULL
      .pmclustEnv$SAVE.iter <- NULL
      .pmclustEnv$CLASS.iter.org <- unlist(apply(.pmclustEnv$Z.spmd, 1,
                                                 which.max))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- try(em.onestep.spmd(PARAM.org))
    if(comm.any(class(PARAM.new) == "try-error")){
      comm.cat("Results of previous iterations are returned.\n", quiet = TRUE)
      .pmclustEnv$CHECK$convergence <- 99
      PARAM.new <- PARAM.org
      break
    }

    .pmclustEnv$CHECK <- check.em.convergence(PARAM.org, PARAM.new, i.iter)
    if(.pmclustEnv$CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      .pmclustEnv$SAVE.param <- c(.pmclustEnv$SAVE.param, PARAM.new)
      CLASS.iter.new <- unlist(apply(.pmclustEnv$Z.spmd, 1, which.max))
      tmp <- as.double(sum(CLASS.iter.new != .pmclustEnv$CLASS.iter.org))
      tmp <- spmd.allreduce.double(tmp, double(1), op = "sum")
      tmp.all <- c(tmp / PARAM.new$N, PARAM.new$logL,
                   PARAM.new$logL - PARAM.org$logL,
                   (PARAM.new$logL - PARAM.org$logL) / PARAM.org$logL)
      .pmclustEnv$SAVE.iter <- rbind(.pmclustEnv$SAVE.iter,
                                     c(tmp, tmp.all, tmp.time))
      .pmclustEnv$CLASS.iter.org <- CLASS.iter.new
    }

    PARAM.org <- PARAM.new
    i.iter <- i.iter + 1
  }

  PARAM.new
} # End of em.step.spmd().

em.onestep.spmd <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "em.Rprof", append = TRUE)
#  }

  PARAM <- m.step.spmd(PARAM)
  e.step.spmd(PARAM)

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step.spmd()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>em.onestep: ", format(Sys.time(), "%H:%M:%S"),
             ", iter: ", .pmclustEnv$CHECK$iter, ", logL: ",
                         sprintf("%-30.15f", PARAM$logL), "\n",
             sep = "", quiet = TRUE)
    if(.pmclustEnv$CONTROL$debug > 4){
      logL <- indep.logL(PARAM)
      comm.cat("  >>indep.logL: ", sprintf("%-30.15f", logL), "\n",
               sep = "", quiet = TRUE)
    }
    if(.pmclustEnv$CONTROL$debug > 20){
      mb.print(PARAM, .pmclustEnv$CHECK)
    }
  }

  PARAM
} # End of em.onestep.spmd().

em.onestep <- em.onestep.spmd

### Obtain classifications.
em.update.class.spmd <- function(){
  .pmclustEnv$CLASS.spmd <- unlist(apply(.pmclustEnv$Z.spmd, 1, which.max))
  invisible()
} # End of em.update.class.spmd().

