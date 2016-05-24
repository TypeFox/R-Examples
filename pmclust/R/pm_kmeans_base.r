### This file provides functions for kmeans.

kmeans.e.step.spmd <- function(PARAM){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  nrow <- nrow(X.spmd)
  ncol <- ncol(X.spmd)

  for(i.k in 1:PARAM$K){
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow, ncol)
    .pmclustEnv$Z.spmd[, i.k] <- sqrt(rowSums(B * B))
  }
  invisible()
} # End of kmeans.e.step.spmd().

kmeans.m.step.spmd <- function(PARAM){
  if(exists("X.spmd", envir = .pmclustEnv)){
    X.spmd <- get("X.spmd", envir = .pmclustEnv)
  }

  for(i.k in 1:PARAM$K){
    id <- .pmclustEnv$CLASS.spmd == i.k
    tmp.n.id <- as.double(sum(id))
    tmp.n.id <- spmd.allreduce.double(tmp.n.id, double(1), op = "sum")

    if(tmp.n.id > 0){
      tmp.sum <- colSums(matrix(X.spmd[id, ], ncol = PARAM$p))
    } else{
      tmp.sum <- rep(0.0, PARAM$p)
    }
    tmp.sum <- spmd.allreduce.double(tmp.sum, double(PARAM$p), op = "sum")

    PARAM$MU[, i.k] <- tmp.sum / tmp.n.id
  } 

  PARAM
} # End of kmeans.m.step.spmd().

kmeans.logL.step.spmd <- function(){
  tmp <- unlist(apply(.pmclustEnv$Z.spmd, 1, which.min))
  tmp.diff <- sum(.pmclustEnv$CLASS.spmd != tmp)

  .pmclustEnv$CLASS.spmd <- tmp
  spmd.allreduce.integer(as.integer(tmp.diff), integer(1), op = "sum")
} # End of kmeans.logL.step.spmd().

check.kmeans.convergence <- function(PARAM.org, PARAM.new, i.iter){
    abs.err <- PARAM.new$logL
    rel.err <- abs.err / PARAM.new$N
    convergence <- 0

    if(i.iter > .pmclustEnv$CONTROL$max.iter){
      convergence <- 2
    } else if(abs.err == 0 || rel.err < .pmclustEnv$CONTROL$rel.err){
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
} # End of check.kmeans.convergence().

kmeans.step.spmd <- function(PARAM.org){
  .pmclustEnv$CHECK <- list(algorithm = "kmeans", i.iter = 0, abs.err = Inf,
                            rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- PARAM.org$N

  ### For debugging.
  if((!is.null(.pmclustEnv$CONTROL$save.log)) && .pmclustEnv$CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .pmclustEnv)){
      .pmclustEnv$SAVE.param <- NULL
      .pmclustEnv$SAVE.iter <- NULL
      .pmclustEnv$CLASS.iter.org <- unlist(apply(.pmclustEnv$Z.spmd, 1,
                                                 which.min))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- kmeans.onestep.spmd(PARAM.org)

    .pmclustEnv$CHECK <- check.kmeans.convergence(PARAM.org, PARAM.new, i.iter)

    if(.pmclustEnv$CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(.pmclustEnv$CONTROL$save.log)) &&
        .pmclustEnv$CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      .pmclustEnv$SAVE.param <- c(.pmclustEnv$SAVE.param, PARAM.new)
      CLASS.iter.new <- unlist(apply(.pmclustEnv$Z.spmd, 1, which.min))
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
} # End of kmeans.step.spmd().

kmeans.onestep.spmd <- function(PARAM){
#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(filename = "kmeans.Rprof", append = TRUE)
#  }

  PARAM <- kmeans.m.step.spmd(PARAM)
  kmeans.e.step.spmd(PARAM)

#  if(.pmclustEnv$COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- kmeans.logL.step.spmd()

  if(.pmclustEnv$CONTROL$debug > 0){
    comm.cat(">>kmeans.onestep: ", format(Sys.time(), "%H:%M:%S"),
             ", iter: ", .pmclustEnv$CHECK$iter, ", logL: ",
                         sprintf("%-30d", PARAM$logL), "\n",
             sep = "", quiet = TRUE)
    if(.pmclustEnv$CONTROL$debug > 10){
      mb.print(PARAM, .pmclustEnv$CHECK)
    }
  }

  PARAM
} # End of kmeans.onestep.spmd().


kmeans.update.class.spmd <- function(){
  .pmclustEnv$CLASS.spmd <- unlist(apply(.pmclustEnv$Z.spmd, 1, which.min))
  invisible()
} # End of kmeans.update.class.spmd().

