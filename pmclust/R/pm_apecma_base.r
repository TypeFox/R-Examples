### This file contains major functions for EM iterations.

### E-step.
apea.step.spmd.k <- function(PARAM, i.k, update.logL = TRUE){
  logdmvnorm(PARAM, i.k)
  update.expectation(PARAM, update.logL = update.logL)
} # End of apea.step.spmd.k().

### CM-step
cm.step.spmd.ETA.MU.SIGMA.k <- function(PARAM, i.k){
  ### MLE For ETA
  PARAM$ETA <- .pmclustEnv$Z.colSums / sum(.pmclustEnv$Z.colSums)
  PARAM$log.ETA <- log(PARAM$ETA)

  ### MLE for MU and SIGMA
  PARAM <- cm.step.spmd.MU.SIGMA.k(PARAM, i.k)

  PARAM
} # End of cm.step.spmd.ETA.MU.SIGMA.k().

cm.step.spmd.MU.SIGMA.k <- function(PARAM, i.k){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  p <- PARAM$p
  p.2 <- p * p

  ### MLE for MU
  tmp.MU <- colSums(X.spmd * .pmclustEnv$Z.spmd[, i.k]) /
            .pmclustEnv$Z.colSums[i.k]
  PARAM$MU[, i.k] <- spmd.allreduce.double(tmp.MU, double(p), op = "sum")

  ### MLE for SIGMA
  if(PARAM$U.check[[i.k]]){
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow(X.spmd), ncol(X.spmd)) *
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

  PARAM
} # End of cm.step.spmd.MU.SIGMA.k().


### APECMa-step.
apecma.step.spmd <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(algorithm = "apecma", i.iter = 0, abs.err = Inf,
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

    PARAM.new <- try(apecma.onestep.spmd(PARAM.org))
    if(comm.any(class(PARAM.new) == "try-error")){
      comm.cat("Results of previous iterations are returned.\n", quiet =TRUE)
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
} # End of apecma.step.spmd().

apecma.onestep.spmd <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "apecma.Rprof", append = TRUE)
#  }

  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.spmd.ETA.MU.SIGMA.k(PARAM, i.k)
    apea.step.spmd.k(PARAM, i.k,
                      update.logL = ifelse(i.k == PARAM$K, TRUE, FALSE))
  }

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- logL.step.spmd()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>apecma.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of apecma.onestep.spmd().

