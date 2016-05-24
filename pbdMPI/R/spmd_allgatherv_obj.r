### Median level functions for R objects. These should not be in S3/S4.

### For general types.
spmd.allgather.object <- function(x,
    comm = .pbd_env$SPMD.CT$comm, unlist = .pbd_env$SPMD.CT$unlist){
  x.raw <- serialize(x, NULL)
  x.count <- spmd.allgather.integer(length(x.raw),
                                    integer(spmd.comm.size(comm = comm)),
                                    comm = comm)
  displs <- c(0L, cumsum(x.count))
  ### mpi_allgatherv() in C requires displs[1:(length(displs) - 1)] only.
  ### It only passes the pointer, so displs[1:length(displs)] is OK, too.
  ret <- spmd.allgatherv.raw(x.raw, raw(sum(x.count)), x.count = x.count,
                             displs = displs, comm = comm)
  ret <- lapply(1:length(x.count),
                function(i) unserialize(ret[(displs[i] + 1):displs[i + 1]]))
  if(unlist){
    ret <- unlist(ret)
  }
  ret
} # End of spmd.allgather.object().

### For array only.
spmd.allgather.array <- function(x,
    comm = .pbd_env$SPMD.CT$comm, unlist = .pbd_env$SPMD.CT$unlist){
  n.dim <- length(dim(x))
  COMM.SIZE <- spmd.comm.size(comm)
  
  all.check <- spmd.allreduce.integer(
                 as.integer(is.double(x) && length(x) > 0),
                 integer(1), op = "sum", comm = comm) == COMM.SIZE
  if(all.check){
    x.dim <- spmd.allgather.integer(dim(x), integer(COMM.SIZE * n.dim),
                                    comm = comm)
    dim(x.dim) <- c(n.dim, COMM.SIZE)
    x.count <- as.integer(apply(x.dim, 2, prod))
    displs <- c(0L, cumsum(x.count))
    ret <- spmd.allgatherv.double(x, double(sum(x.count)), x.count,
                                  displs = displs, comm = comm)
    ret <- lapply(1:COMM.SIZE,
                  function(i){
                    array(data = ret[(displs[i] + 1):displs[i + 1]],
                          dim = x.dim[, i])
                  })
    if(unlist){
      ret <- unlist(ret)
    }
    return(ret)
  }

  all.check <- spmd.allreduce.integer(
                 as.integer(is.integer(x) && length(x) > 0),
                 integer(1), op = "sum", comm = comm) == COMM.SIZE
  if(all.check){
    x.dim <- spmd.allgather.integer(dim(x), integer(COMM.SIZE * n.dim),
                                    comm = comm)
    dim(x.dim) <- c(n.dim, COMM.SIZE)
    x.count <- as.integer(apply(x.dim, 2, prod))
    displs <- c(0L, cumsum(x.count))
    ret <- spmd.allgatherv.integer(x, integer(sum(x.count)), x.count,
                                   displs = displs, comm = comm)
    ret <- lapply(1:COMM.SIZE,
                  function(i){
                    array(data = ret[(displs[i] + 1):displs[i + 1]],
                          dim = x.dim[, i])
                  })
    if(unlist){
      ret <- unlist(ret)
    }
    return(ret)
  }

  spmd.allgather.object(x, comm = comm, unlist = unlist)
} # End of spmd.allgather.array().

