### This file contains functions to load balance of data X.gbd.
### Assume gbd.major = 1.

comm.balance.info <- function(X.gbd,
    balance.method = .pbd_env$SPMD.IO$balance.method[1],
    comm = .pbd_env$SPMD.CT$comm){
  ### Check gbd.
  if(!comm.allcommon.integer(length(dim(X.gbd)), comm = comm)){
    comm.stop("Dimension of X.gbd should all equal to 2.", comm = comm)
  }
  if(!comm.allcommon.integer(ncol(X.gbd), comm = comm)){
    comm.stop("X.gbd should have the same # of columns.", comm = comm)
  }

  COMM.SIZE <- spmd.comm.size(comm)
  COMM.RANK <- spmd.comm.rank(comm)

  N.gbd <- nrow(X.gbd)
  N.allgbd <- spmd.allgather.integer(as.integer(N.gbd), integer(COMM.SIZE),
                                     comm = comm)
  N <- sum(N.allgbd)

  ### Build table by method.
  if(balance.method[1] == "block"){
    n <- floor(N / COMM.SIZE)
    n.residual <- N %% COMM.SIZE
    new.N.allgbd <- rep(n, COMM.SIZE) +
                    rep(c(0, 1), c(COMM.SIZE - n.residual, n.residual))
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allgbd)
  } else if(balance.method[1] == "block0"){
    n <- floor(N / COMM.SIZE)
    n.residual <- N %% COMM.SIZE
    new.N.allgbd <- rep(n, COMM.SIZE) +
                    rep(c(1, 0), c(n.residual, COMM.SIZE - n.residual))
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allgbd)
  } else if(balance.method[1] == "block.cyclic"){
    n <- ceiling(N / COMM.SIZE)
    rep.n <- N %/% n
    new.N.allgbd <- rep(n, rep.n)
    if(n * rep.n < N){
      new.N.allgbd <- c(new.N.allgbd, (N - n * rep.n))
    }
    if(length(new.N.allgbd) < COMM.SIZE){
      new.N.allgbd <- c(new.N.allgbd,
                        rep(0, COMM.SIZE - length(new.N.allgbd)))
    }
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allgbd)
  } else{
    comm.stop("balance.method is undefined.", comm = comm)
  }

  rank.org <- rep(0:(COMM.SIZE - 1), N.allgbd)

  ### Build send and recv information if any.
  send.info <- data.frame(org = rank.org[rank.org == COMM.RANK],
                          belong = rank.belong[rank.org == COMM.RANK])
  recv.info <- data.frame(org = rank.org[rank.belong == COMM.RANK],
                          belong = rank.belong[rank.belong == COMM.RANK])

  list(send = send.info, recv = recv.info, N.allgbd = N.allgbd,
       new.N.allgbd = new.N.allgbd, balance.method = balance.method[1])
} # End of comm.balance.info().


comm.load.balance <- function(X.gbd, bal.info = NULL,
    balance.method = .pbd_env$SPMD.IO$balance.method[1],
    comm = .pbd_env$SPMD.CT$comm){
  ### Check.
  if(!comm.allcommon.integer(length(dim(X.gbd)), comm = comm)){
    comm.stop("Dimension of X.gbd should all equal to 2.", comm = comm)
  }
  if(!comm.allcommon.integer(ncol(X.gbd), comm = comm)){
    comm.stop("X.gbd should have the same # of columns.", comm = comm)
  }

  ### Get bal.info if NULL.
  COMM.RANK <- spmd.comm.rank(comm)
  if(is.null(bal.info)){
    bal.info <- comm.balance.info(X.gbd, balance.method, comm = comm)
  }

  p <- ncol(X.gbd)

  ### Redistributing.
  send.to <- as.integer(unique(bal.info$send$belong))
  if(length(send.to) > 0){
    for(i in send.to){
      if(i != COMM.RANK){
        tmp <- X.gbd[bal.info$send$belong == i,]
        spmd.isend.default(tmp, rank.dest = i, tag = COMM.RANK, comm = comm)
      }
    }
  }

  recv.from <- as.integer(unique(bal.info$recv$org))
  if(length(recv.from) > 0){
    ret <- NULL
    for(i in recv.from){
      if(i != COMM.RANK){
        tmp <- spmd.recv.default(rank.source = i, tag = i, comm = comm)
      } else{
        tmp <- X.gbd[bal.info$send$belong == i,]
      }
      colnames(tmp) <- colnames(X.gbd)
      ret <- base::rbind(ret, tmp)
    }
  } else{
    ret <- X.gbd
  }

  if(bal.info$new.N.allgbd[spmd.comm.rank(comm) + 1] == 0){
    ret <- ret[0,]
  }

  spmd.wait()

  ### Return.
  ret
} # End of comm.load.balance().


comm.unload.balance <- function(new.X.gbd, bal.info,
    comm = .pbd_env$SPMD.CT$comm){
  rev.bal.info <- list(send = data.frame(org = bal.info$recv$belong,
                                         belong = bal.info$recv$org),
                       recv = data.frame(org = bal.info$send$belong,
                                         belong = bal.info$send$org),
                       N.allgbd = bal.info$new.N.allgbd,
                       new.N.allgbd = bal.info$N.allgbd)
  comm.load.balance(new.X.gbd, bal.info = rev.bal.info, comm = comm)
} # End of comm.unload.balance().

