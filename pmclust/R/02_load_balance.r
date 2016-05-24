### This file contains functions to load balance of data X.spmd.

balance.info <- function(X.spmd, comm = .pbd_env$SPMD.CT$comm,
    spmd.major = 1, method = c("block.cyclic", "block0")){
  COMM.SIZE <- spmd.comm.size(comm)
  COMM.RANK <- spmd.comm.rank(comm)

  if(!is.matrix(X.spmd)){
    X.spmd <- as.matrix(X.spmd)
  }

  if(spmd.major == 1){
    N.spmd <- nrow(X.spmd)
  } else if(spmd.major == 2){
    N.spmd <- ncol(X.spmd)
  } else{
    stop("spmd.major = 1 or 2.")
  }
  N.allspmd <- spmd.allgather.integer(as.integer(N.spmd), integer(COMM.SIZE),
                                      comm = comm)
  N <- sum(N.allspmd)

  if(method[1] == "block.cyclic"){
    ### This can cause problems.
    # n <- ceiling(N / COMM.SIZE)
    # new.N.allspmd <- c(rep(n, COMM.SIZE - 1), N - n * (COMM.SIZE - 1))
    # rank.org <- rep(0:(COMM.SIZE - 1), N.allspmd)
    # rank.belong <- rep(0:(COMM.SIZE - 1), each = n)[1:N]

    ### Try again.
    n <- ceiling(N / COMM.SIZE)
    rep.n <- N %/% n
    new.N.allspmd <- rep(n, rep.n)
    if(n * rep.n < N){
      new.N.allspmd <- c(new.N.allspmd, (N - n * rep.n))
    }
    if(length(new.N.allspmd) < COMM.SIZE){
      new.N.allspmd <- c(new.N.allspmd,
                         rep(0, COMM.SIZE - length(new.N.allspmd)))
    }
    rank.org <- rep(0:(COMM.SIZE - 1), N.allspmd)
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allspmd) 
  } else if(method[1] == "block0"){
    ### Try block0 method which is a better way to balance data. However,
    ### this is not necessary in block-cyclic, so useless for ddmatrix.
    n <- floor(N / COMM.SIZE)
    n.residual <- N %% COMM.SIZE
    new.N.allspmd <- rep(n, COMM.SIZE) +
                     rep(c(1, 0), c(n.residual, COMM.SIZE - n.residual))
    rank.org <- rep(0:(COMM.SIZE - 1), N.allspmd)
    rank.belong <- rep(0:(COMM.SIZE - 1), new.N.allspmd)
  } else{
    comm.stop("method is not found.")
  }

  ### Build send and recv information if any.
  send.info <- data.frame(org = rank.org[rank.org == COMM.RANK],
                          belong = rank.belong[rank.org == COMM.RANK])
  recv.info <- data.frame(org = rank.org[rank.belong == COMM.RANK],
                          belong = rank.belong[rank.belong == COMM.RANK])

  list(send = send.info, recv = recv.info, N.allspmd = N.allspmd,
       new.N.allspmd = new.N.allspmd, spmd.major = spmd.major)
} # End of balance.info()


load.balance <- function(X.spmd, bal.info = NULL, comm = .pbd_env$SPMD.CT$comm,
    spmd.major = 1){
  COMM.RANK <- spmd.comm.rank(comm)
  if(is.null(bal.info)){
    bal.info <- balance.info(X.spmd, comm = comm, spmd.major = spmd.major)
  }

  if(!is.matrix(X.spmd)){
    X.spmd <- as.matrix(X.spmd)
  }
  if(spmd.major == 1){
    p <- ncol(X.spmd)
  } else if(spmd.major == 2){
    p <- nrow(X.spmd)
  } else{
    stop("spmd.major = 1 or 2.")
  }

  storage.mode(X.spmd) <- "double"

  send.to <- as.integer(unique(bal.info$send$belong))
  if(length(send.to) > 0){
    if(spmd.major == 1){
      for(i in send.to){
        if(i != COMM.RANK){
          tmp <- matrix(X.spmd[bal.info$send$belong == i,], ncol = p)
          spmd.isend.double(tmp, rank.dest = i, tag = COMM.RANK, comm = comm)
        }
      }
    } else{
      for(i in send.to){
        if(i != COMM.RANK){
          tmp <- matrix(X.spmd[, bal.info$send$belong == i], nrow = p)
          spmd.isend.double(tmp, rank.dest = i, tag = COMM.RANK, comm = comm)
        }
      }
    }
  }

  recv.from <- as.integer(unique(bal.info$recv$org))
  if(length(recv.from) > 0){
    ret <- NULL
    if(spmd.major == 1){
      for(i in recv.from){
        if(i != COMM.RANK){
          total.row <- sum(bal.info$recv$org == i)
          tmp <- spmd.recv.double(double(total.row * p),
                                  rank.source = i, tag = i, comm = comm)
          dim(tmp) <- c(total.row, p)
        } else{
          tmp <- matrix(X.spmd[bal.info$send$belong == i,], ncol = p)
        }
        ret <- base::rbind(ret, tmp)
      }
    } else{
      for(i in recv.from){
        if(i != COMM.RANK){
          total.column <- sum(bal.info$recv$org == i)
          tmp <- spmd.recv.double(double(total.column * p),
                                  rank.source = i, tag = i, comm = comm)
          dim(tmp) <- c(p, total.column)
        } else{
          tmp <- matrix(X.spmd[, bal.info$send$belong == i], nrow = p)
        }
        ret <- base::cbind(ret, tmp)
      }
    }
  } else{
    ret <- X.spmd
  }

  if(bal.info$new.N.allspmd[spmd.comm.rank(comm) + 1] == 0){
    if(spmd.major == 1){
      ret <- matrix(0, nrow = 0, ncol = p)
    } else{
      ret <- matrix(0, nrow = p, ncol = 0)
    }
  }

  spmd.wait()

  ret
} # End of load.balance().


unload.balance <- function(new.X.spmd, bal.info, comm = .pbd_env$SPMD.CT$comm){
  rev.bal.info <- list(send = data.frame(org = bal.info$recv$belong,
                                         belong = bal.info$recv$org),
                       recv = data.frame(org = bal.info$send$belong,
                                         belong = bal.info$send$org),
                       N.allspmd = bal.info$new.N.allspmd,
                       new.N.allspmd = bal.info$N.allspmd,
                       spmd.major = bal.info$spmd.major)
  load.balance(new.X.spmd, bal.info = rev.bal.info, comm = comm)
} # End of unload.balance().

