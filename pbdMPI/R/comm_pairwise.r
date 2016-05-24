### Distributed all pair-wise functions.
### Assume gbd.major = 1.

comm.pairwise <- function(X, pairid.gbd = NULL,
    FUN = function(x, y, ...){ return(as.vector(dist(rbind(x, y), ...))) },
    ..., diag = FALSE, symmetric = TRUE, comm = .pbd_env$SPMD.CT$comm){
  ### Check if all non-NULL.
  check <- spmd.allreduce.integer(!is.null(pairid.gbd), integer(1),
                                  comm = comm) == spmd.comm.size(comm)

  ### Evaluate accordingly.
  if(check){
    ret <- comm.pairwise.common(X, pairid.gbd, FUN, ..., diag = diag,
                                symmetric = symmetric, comm = comm)
  } else{
    ret <- comm.pairwise.gbd(X, FUN, ..., diag = diag,
                             symmetric = symmetric, comm = comm)
  }

  ret
} # End of comm.pairwise().

### Input a common matrix and a gbd pairid.
comm.pairwise.common <- function(X, pairid.gbd,
    FUN = function(x, y, ...){ return(as.vector(dist(rbind(x, y), ...))) },
    ..., diag = FALSE, symmetric = TRUE, comm = .pbd_env$SPMD.CT$comm){
  ### FUN <- function(x, y, ...) is a user defined function.

  ### Check pairid.gbd.
  if(!comm.allcommon.integer(length(dim(pairid.gbd)), comm = comm)){
    comm.stop("Dimension of pairid.gbd should all equal to 2.", comm = comm)
  }
  if(!comm.allcommon.integer(ncol(pairid.gbd), comm = comm)){
    comm.stop("pairid.gbd should have the same # of columns.", comm = comm)
  }

  ### Evaluate the user defined function
  if(nrow(pairid.gbd) > 0){
    ret <- lapply(1:nrow(pairid.gbd),
             function(i.pair){
               FUN(X[pairid.gbd[i.pair, 1],], X[pairid.gbd[i.pair, 2],], ...)
             }, ...) 
    ret <- cbind(pairid.gbd, do.call("c", ret))
  } else{
    ret <- matrix(0, nrow = 0, ncol = 3)
  }

  ### Return.
  colnames(ret) <- c("i", "j", "value")
  rownames(ret) <- NULL
  ret
} # End of comm.pairwise.common().

### Input a gbd matrix.
comm.pairwise.gbd <- function(X.gbd,
    FUN = function(x, y, ...){ return(as.vector(dist(rbind(x, y), ...))) },
    ..., diag = FALSE, symmetric = TRUE, comm = .pbd_env$SPMD.CT$comm){
  ### FUN <- function(x, y, ...) is a user defined function.

  ### Check X.gbd.
  if(!comm.allcommon.integer(length(dim(X.gbd)), comm = comm)){
    comm.stop("Dimension of X.gbd should all equal to 2.", comm = comm)
  }
  if(!comm.allcommon.integer(ncol(X.gbd), comm = comm)){
    comm.stop("X.gbd should have the same # of columns.", comm = comm)
  }

  ### Get info.
  COMM.RANK <- spmd.comm.rank(comm)
  COMM.SIZE <- spmd.comm.size(comm)

  N.gbd <- nrow(X.gbd)
  N.allgbd <- spmd.allgather.integer(as.integer(N.gbd), integer(COMM.SIZE),
                                     comm = comm)
  N <- sum(N.allgbd)
  N.cumsum <- c(1, cumsum(N.allgbd) + 1)

  ### Only local diagonal block.
  ret.local <- matrix(0.0, nrow = N.gbd, ncol = N.gbd)
  if(N.gbd > 0){
    for(i in 1:N.gbd){
      for(j in 1:N.gbd){
        ### Check.
        flag <- FALSE
        if(i > j){                      ### lower-triangular.
          flag <- TRUE
        }
        if((!symmetric) && (i < j)){    ### Upper-triangular.
          flag <- TRUE
        }
        if(diag && (i == j)){           ### Diagonals.
          flag <- TRUE
        }
        if(flag){                       ### Compute.
          ret.local[i, j] <- FUN(X.gbd[i,], X.gbd[j,], ...)
        }
      }
    }
  }

  #### Obtain all other ranks if more than one processor.
  if(COMM.SIZE > 1){
    ret.lower <- NULL
    ret.upper <- NULL

    ### Evaluate the user defined function
    for(i.rank in 0:(COMM.SIZE - 1)){
      if(N.allgbd[i.rank + 1] != 0){
        X.other <- bcast(X.gbd, rank.source = i.rank, comm = comm)

        if(N.gbd > 0){
          ### Lower-triangular block.
          if(COMM.RANK < i.rank){
            tmp <- matrix(0.0, nrow = N.allgbd[i.rank + 1], ncol = N.gbd)
            for(i in 1:N.allgbd[i.rank + 1]){
              for(j in 1:N.gbd){
                tmp[i, j] <- FUN(X.other[i,], X.gbd[j,], ...)
              }
            }
            ret.lower <- rbind(ret.lower, tmp)
          }

          ### Upper-triangular block.
          if(COMM.RANK > i.rank && !symmetric){
            tmp <- matrix(0.0, nrow = N.allgbd[i.rank + 1], ncol = N.gbd)
            for(i in 1:N.allgbd[i.rank + 1]){
              for(j in 1:N.gbd){
                tmp[i, j] <- FUN(X.other[i,], X.gbd[j,], ...)
              }
            }
            ret.upper <- rbind(ret.upper, tmp)
          }
        }
      }
    }

    ### Combine all blocks.
    if(N.gbd > 0){
      ret.local <- rbind(ret.upper, ret.local, ret.lower)
    }
  }
  dim(ret.local) <- c(length(ret.local), 1)

  ### Add and combine i, j, and value.
  if(N.gbd > 0){
    if(symmetric){
      ret <- cbind(rep(N.cumsum[COMM.RANK + 1]:N, N.gbd),
                   rep(N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1),
                       each = N - N.cumsum[COMM.RANK + 1] + 1),
                   ret.local)
      ret <- ret[ret[, 1] >= ret[, 2],]
    } else{
      ret <- cbind(rep(1:N, N.gbd),
                   rep(N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1),
                       each = N),
                   ret.local)
    }
    dim(ret) <- c(length(ret) / 3, 3)       ### Could be all gone.

    if(!diag && nrow(ret) > 0){
      ret <- ret[ret[, 1] != ret[, 2],]
    }

    dim(ret) <- c(length(ret) / 3, 3)       ### Could be all gone.
  } else{
    ret <- matrix(0, nrow = 0, ncol = 3)    ### Could be all gone.
  }

  ### Return.
  colnames(ret) <- c("i", "j", "value")
  rownames(ret) <- NULL
  ret
} # End of comm.pairwise.gbd().

