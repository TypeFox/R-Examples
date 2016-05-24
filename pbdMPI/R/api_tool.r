### This file contains some useful tools.

divide.job.vector <- function(n, method, COMM.SIZE, COMM.RANK){
  jid <- NULL
  if(n > COMM.SIZE){
    if(method[1] == "block"){
      n.job <- floor(n / COMM.SIZE)
      n.residual <- n %% COMM.SIZE
      if(n.residual == 0){
        jid <- 1:n.job + COMM.RANK * n.job
      } else{
        if(COMM.RANK < (COMM.SIZE - n.residual)){
          jid <- 1:n.job + COMM.RANK * n.job
        } else{
          jid <- 1:(n.job + 1) + COMM.RANK * (n.job + 1) -
                 (COMM.SIZE - n.residual)
        }
      }
    } else if(method[1] == "block0"){
      n.job <- floor(n / COMM.SIZE)
      n.residual <- n %% COMM.SIZE
      if(COMM.RANK < n.residual){
        jid <- 1:(n.job + 1) + COMM.RANK * (n.job + 1)
      } else{
        jid <- 1:n.job + n.residual * (n.job + 1) +
               (COMM.RANK - n.residual) * n.job
      }
    } else if(method[1] == "cycle"){
      jid <- which((0:(n - 1)) %% COMM.SIZE == COMM.RANK)
    }
  } else{
    if(COMM.RANK < n){
      jid <- COMM.RANK + 1
    }
  }

  jid
} # End of divide.job.vector().

divide.job.list <- function(n, method, COMM.SIZE, COMM.RANK){
  jid <- list(start = NULL, end = NULL, by = 0, length = 0)
  if(n > COMM.SIZE){
    if(method[1] == "block"){
      n.job <- floor(n / COMM.SIZE)
      n.residual <- n %% COMM.SIZE

      if(n.residual == 0){
        jid <- list(start = 1 + COMM.RANK * n.job,
                    end = n.job + COMM.RANK * n.job,
                    by = 1,
                    length = n.job)
      } else{
        if(COMM.RANK < (COMM.SIZE - n.residual)){
          jid <- list(start = 1 + COMM.RANK * n.job,
                      end = n.job + COMM.RANK * n.job,
                      by = 1,
                      length = n.job)
        } else{
          jid <- list(start = 1 + COMM.RANK * (n.job + 1) -
                              (COMM.SIZE - n.residual),
                      end = (n.job + 1) + COMM.RANK * (n.job + 1) -
                            (COMM.SIZE - n.residual),
                      by = 1,
                      length = n.job + 1)
        }
      }
    } else if(method[1] == "block0"){
      n.job <- floor(n / COMM.SIZE)
      n.residual <- n %% COMM.SIZE
      if(COMM.RANK < n.residual){
        jid <- list(start = 1 + COMM.RANK * (n.job + 1),
                    end = (n.job + 1) + COMM.RANK * (n.job + 1),
                    by = 1,
                    length = n.job + 1)
      } else{
        jid <- list(start = 1 + n.residual * (n.job + 1) +
                            (COMM.RANK - n.residual) * n.job,
                    end = n.job + n.residual * (n.job + 1) +
                          (COMM.RANK - n.residual) * n.job,
                    by = 1,
                    length = n.job)
      }
    } else if(method[1] == "cycle"){
      tl.cycle <- n %/% COMM.SIZE
      tl.remainder <- n %% COMM.SIZE
      my.remainder <- (COMM.RANK < tl.remainder)
      jid <- list(start = COMM.RANK + 1,
                  end = (tl.cycle - 1 + my.remainder) * COMM.SIZE +
                        COMM.RANK + 1,
                  by = COMM.SIZE,
                  length = tl.cycle + my.remainder)
    }
  } else{
    if(COMM.RANK < n){
      jid <- list(start = COMM.RANK + 1,
                  end = COMM.RANK + 1,
                  by = 0,
                  length = 1)
    }
  }

  ### Force by = 0 if only one element.
  if(jid$length == 1){
    jid$by <- 0
  }

  jid
} # End of divide.job.list().

divide.job <- function(n, method = .pbd_env$SPMD.CT$divide.method[1],
    comm = .pbd_env$SPMD.CT$comm, reduced = FALSE){
  COMM.SIZE <- spmd.comm.size(comm)
  COMM.RANK <- spmd.comm.rank(comm)

  if(!reduced){
    jid <- divide.job.vector(n, method, COMM.SIZE, COMM.RANK)
  } else{
    jid <- divide.job.list(n, method, COMM.SIZE, COMM.RANK)
  }

  jid
} # End of divide.job().

divide.job.all <- function(n, method = .pbd_env$SPMD.CT$divide.method[1],
    comm = .pbd_env$SPMD.CT$comm, reduced = FALSE){
  COMM.SIZE <- spmd.comm.size(comm)

  alljid <- rep(list(NULL), COMM.SIZE)
  if(!reduced){
    for(i.rank in 1:COMM.SIZE){
      alljid[[i.rank]] <- divide.job.vector(n, method, COMM.SIZE, i.rank - 1)
    }
  } else{
    for(i.rank in 1:COMM.SIZE){
      alljid[[i.rank]] <- divide.job.list(n, method, COMM.SIZE, i.rank - 1)
    }
  }

  alljid
} # End of divide.job.all().

get.jid <- function(n, method = .pbd_env$SPMD.CT$divide.method[1],
    all = FALSE,
    comm = .pbd_env$SPMD.CT$comm, reduced = FALSE){
  if(! method[1] %in% .pbd_env$SPMD.CT$divide.method){
    stop("The method for dividing jobs is not found.")
  }

  if(all){
    ret <- divide.job.all(n, method = method, comm = comm, reduced = reduced)
  } else{
    ret <- divide.job(n, method = method, comm = comm, reduced = reduced)
  }

  ret
} # End of get.jid().

