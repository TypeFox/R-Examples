### This file contains major functions for EM iterations.

### E-step.
ape.step.spmd <- function(PARAM){
  for(i.k in 1:PARAM$K){
    logdmvnorm(PARAM, i.k)
  }

  ape.update.expectation(PARAM)
} # End of ape.step.spmd().

ape.step.spmd.k <- function(PARAM, i.k, update.logL = TRUE){
  logdmvnorm(PARAM, i.k)
  ape.update.expectation.k(PARAM, i.k, update.logL)
} # End of ape.step.spmd.k().


### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
ape.update.expectation <- function(PARAM, update.logL = TRUE){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  N <- nrow(X.spmd)
  K <- PARAM$K

  .pmclustEnv$W.spmd <- W.plus.y(.pmclustEnv$W.spmd, PARAM$log.ETA, N, K)
  .pmclustEnv$U.spmd <- exp(.pmclustEnv$W.spmd)
  .pmclustEnv$Z.spmd <- .pmclustEnv$U.spmd

  tmp.id <- rowSums(.pmclustEnv$W.spmd < .pmclustEnv$CONTROL$exp.min) == K |
            rowSums(.pmclustEnv$W.spmd > .pmclustEnv$CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.spmd <- .pmclustEnv$W.spmd[tmp.id,]

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
} # End of ape.update.expectation().

ape.update.expectation.k <- function(PARAM, i.k, update.logL = TRUE){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  N <- nrow(X.spmd)
  K <- PARAM$K

  .pmclustEnv$W.spmd[, i.k] <- W.plus.y.k(.pmclustEnv$W.spmd, PARAM$log.ETA,
                                          N, K, i.k)
  .pmclustEnv$U.spmd[, i.k] <- exp(.pmclustEnv$W.spmd[, i.k])
  .pmclustEnv$Z.spmd <- .pmclustEnv$U.spmd

  tmp.id <- rowSums(.pmclustEnv$W.spmd < .pmclustEnv$CONTROL$exp.min) == K |
            rowSums(.pmclustEnv$W.spmd > .pmclustEnv$CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.spmd <- .pmclustEnv$W.spmd[tmp.id,]

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
} # End of ape.update.expectation.k().


### APECM-step.
apecm.step.spmd <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(algorithm = "apecm", i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
      .pmclustEnv$CONTROL$save.log){
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

    PARAM.new <- try(apecm.onestep.spmd(PARAM.org))
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
} # End of apecm.step.spmd().

apecm.onestep.spmd <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "apecm.Rprof", append = TRUE)
#  }

  ### Update ETA
  PARAM <- cm.step.spmd.ETA(PARAM)
  ape.step.spmd(PARAM)

  ### Update MU and SIGMA
  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.spmd.MU.SIGMA.k(PARAM, i.k)
    ape.step.spmd.k(PARAM, i.k,
                       update.logL = ifelse(i.k == PARAM$K, TRUE, FALSE))
  }

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step.spmd()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>apecm.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of apecm.onestep.spmd().

