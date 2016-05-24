### Distributed distance.
### Assume gbd.major = 1.

comm.dist <- function(X.gbd, method = "euclidean", diag = FALSE,
    upper = FALSE, p = 2, comm = .pbd_env$SPMD.CT$comm,
    return.type = c("common", "gbd")){
  if(return.type[1] == "common"){
    ret <- comm.dist.common(X.gbd, method = method, diag = diag, upper = upper,
                            p = p, comm = comm)
  } else if(return.type[1] == "gbd"){
    ret <- comm.dist.gbd(X.gbd, method = method, diag = diag, upper = upper,
                         p = p, comm = comm)
  } else{
    comm.stop("The return.type is not found.", comm = comm)
  }
  ret
} # End of comm.dist().

### Return a common distance matrix.
comm.dist.common <- function(X.gbd, method = "euclidean", diag = FALSE,
    upper = FALSE, p = 2, comm = .pbd_env$SPMD.CT$comm){
  ### Check.
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

  ### Allocate a full distance matrix.
  ret <- matrix(-Inf, nrow = N, ncol = N)

  ### Only local diagonal block.
  tmp <- as.matrix(dist(X.gbd, method = method,
                        diag = diag, upper = upper, p = p))
  ret[N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1),
      N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1)] <- tmp

  ### Obtain all other ranks distance.
  for(i.rank in 0:(COMM.SIZE - 1)){
    if(N.allgbd[i.rank + 1] != 0){
      X.other <- bcast(X.gbd, rank.source = i.rank, comm = comm)

      if(COMM.RANK < i.rank){
        ### Keep the right order.
        tmp <- as.matrix(dist(rbind(X.gbd, X.other), method = method,
                              diag = diag, upper = upper, p = p))

        ### Replace the lower-triangular block.
        ret[N.cumsum[i.rank + 1]:(N.cumsum[i.rank + 2] - 1),
            N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1)] <-
          tmp[(nrow(X.gbd) + 1):nrow(tmp), 1:nrow(X.gbd)]
      } else if(COMM.RANK > i.rank){
        ### Keep the right order.
        tmp <- as.matrix(dist(rbind(X.other, X.gbd), method = method,
                              diag = diag, upper = upper, p = p))

        ### Replace the lower-triangular block.
        ret[N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1),
            N.cumsum[i.rank + 1]:(N.cumsum[i.rank + 2] - 1)] <-
          tmp[(nrow(X.other) + 1):nrow(tmp), 1:nrow(X.other)]
      }
    }
  }

  ### Return.
  ret <- allreduce(ret, op = "max")
  ret <- as.dist(ret, diag = diag, upper = upper)
  ret
} # End of comm.dist.common().

### Return a gbd distance matrix.
comm.dist.gbd <- function(X.gbd, method = "euclidean", diag = FALSE,
    upper = FALSE, p = 2, comm = .pbd_env$SPMD.CT$comm){
  ### Check.
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
  if(N.gbd > 0){
    ret <- as.matrix(dist(X.gbd, method = method, diag = diag, upper = upper,
                          p = p))
  }

  ### Obtain all other ranks distance.
  if(COMM.SIZE > 1){
    for(i.rank in 1:(COMM.SIZE - 1)){
      if(N.allgbd[i.rank + 1] != 0){
        X.other <- bcast(X.gbd, rank.source = i.rank, comm = comm)

        if(COMM.RANK < i.rank && N.gbd > 0){
          ### Keep the right order.
          tmp <- as.matrix(dist(rbind(X.gbd, X.other), method = method,
                                diag = diag, upper = upper, p = p))

          ### Only need the lower-triangular block.
          ret <- rbind(ret, tmp[(nrow(X.gbd) + 1):nrow(tmp), 1:nrow(X.gbd)])
        }
      }
    }
  }

  ### Summarized for returning.
  if(N.gbd > 0){
    ### Build indices
    ret <- cbind(rep(N.cumsum[COMM.RANK + 1]:N, N.allgbd[COMM.RANK + 1]),
                 rep(N.cumsum[COMM.RANK + 1]:(N.cumsum[COMM.RANK + 2] - 1),
                     each = N - N.cumsum[COMM.RANK + 1] + 1),
                 as.vector(ret))

    ### Trim the diagonals and upper triangles.
    ret <- ret[ret[, 1] > ret[, 2],]
  } else{
    ret <- matrix(0, nrow = 0, ncol = 3)
  }

  ### Return.
  dim(ret) <- c(length(ret) / 3, 3)
  colnames(ret) <- c("i", "j", "value")
  rownames(ret) <- NULL
  ret
} # End of comm.dist.gbd().
