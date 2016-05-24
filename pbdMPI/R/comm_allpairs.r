### Distributed all pairs functions.
### Assume gbd.major = 1.

comm.allpairs <- function(N, diag = FALSE, symmetric = TRUE,
    comm = .pbd_env$SPMD.CT$comm){
  ### Check.
  if(!comm.allcommon.integer(N, comm = comm)){
    comm.stop("N should be all the same.")
  }

  ### Generate index.matrix.
  jid <- get.jid(N * N)
  ret <- cbind(rep(1:N, N)[jid], rep(1:N, each = N)[jid])
  if(length(ret) > 0 && symmetric){
    dim(ret) <- c(length(ret) / 2, 2)
    ret <- ret[ret[, 1] >= ret[, 2],]
  }
  if(length(ret) > 0 && !diag){
    dim(ret) <- c(length(ret) / 2, 2)
    ret <- ret[ret[, 1] != ret[, 2],]
  }

  ### Check.
  if(length(ret) > 0){
    dim(ret) <- c(length(ret) / 2, 2)
  } else{
    ret <- matrix(0, nrow = 0, ncol = 2)
  }

  ### Return.
  ret <- comm.load.balance(ret, comm = comm)
  colnames(ret) <- c("i", "j")
  rownames(ret) <- NULL
  ret
} # End of comm.allpairs().
