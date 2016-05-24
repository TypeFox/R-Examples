### Median level functions for R objects. These should not be in S3/S4.

### For general types.
spmd.scatter.object <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  COMM.RANK <- spmd.comm.rank(comm)

  if(COMM.RANK == rank.source){
    x.raw <- lapply(1:spmd.comm.size(comm), function(i) serialize(x[[i]], NULL))
    x.count <- spmd.bcast.integer(do.call("c", lapply(x.raw, length)),
                                  rank.source = rank.source, comm = comm)
    ret <- unserialize(spmd.scatterv.raw(do.call("c", x.raw),
                                         raw(x.count[COMM.RANK + 1]),
                                         x.count = x.count,
                                         rank.source = rank.source,
                                         comm = comm))
  } else{
    x.count <- spmd.bcast.integer(integer(spmd.comm.size(comm)),
                                  rank.source = rank.source, comm = comm)
    ret <- unserialize(spmd.scatterv.raw(raw(0),
                                         raw(x.count[COMM.RANK + 1]),
                                         x.count = x.count,
                                         rank.source = rank.source,
                                         comm = comm))
  }
  ret
} # End of spmd.scatter.object().

### For array only.
spmd.scatter.array <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- spmd.comm.size(comm)
  COMM.RANK <- spmd.comm.rank(comm)

  if(COMM.RANK == rank.source){
    check.double <- do.call("sum",
                      lapply(1:COMM.SIZE,
                             function(i) is.double(x[[i]]))) == COMM.SIZE
    check.integer <- do.call("sum",
                       lapply(1:COMM.SIZE,
                              function(i) is.integer(x[[i]]))) == COMM.SIZE
    check <- c(check.double, check.integer)
    spmd.bcast.integer(as.integer(check), rank.source = rank.source,
                       comm = comm)
  } else{
    check <- spmd.bcast.integer(integer(2), rank.source = rank.source,
                                comm = comm)
    check.double <- check[1]
    check.integer <- check[2]
  }

  if(check.double || check.integer){
    if(COMM.RANK == rank.source){
      x.dim <- lapply(x, dim)
      x.dim.count <- spmd.bcast.integer(do.call("c", lapply(x.dim, length)),
                                        rank.source = rank.source, comm = comm)
      ret.dim <- spmd.scatterv.integer(do.call("c", x.dim),
                                       integer(x.dim.count[COMM.RANK + 1]),
                                       x.count = x.dim.count,
                                       rank.source = rank.source,
                                       comm = comm)
    } else{
      x.dim.count <- spmd.bcast.integer(integer(COMM.SIZE),
                                        rank.source = rank.source, comm = comm)
      ret.dim <- spmd.scatterv.integer(integer(0),
                                       integer(x.dim.count[COMM.RANK + 1]),
                                       x.count = x.dim.count,
                                       rank.source = rank.source,
                                       comm = comm)
    }
  } else{
    ret <- spmd.scatter.object(x, rank.source = rank.source, comm = comm)
    return(ret)
  }

  if(check.double){
    if(COMM.RANK == rank.source){
      x.count <- spmd.bcast.integer(do.call("c", lapply(x, length)),
                                    rank.source = rank.source, comm = comm)
      ret <- spmd.scatterv.double(do.call("c", x),
                                  double(x.count[COMM.RANK + 1]),
                                  x.count = x.count,
                                  rank.source = rank.source,
                                  comm = comm)
    } else{
      x.count <- spmd.bcast.integer(integer(COMM.SIZE),
                                    rank.source = rank.source, comm = comm)
      ret <- spmd.scatterv.double(double(0),
                                  double(x.count[COMM.RANK + 1]),
                                  x.count = x.count,
                                  rank.source = rank.source,
                                  comm = comm)
    }
    dim(ret) <- ret.dim
    return(ret)
  }

  if(check.integer){
    if(COMM.RANK == rank.source){
      x.count <- spmd.bcast.integer(do.call("c", lapply(x, length)),
                                    rank.source = rank.source, comm = comm)
      ret <- spmd.scatterv.integer(do.call("c", x),
                                   integer(x.count[COMM.RANK + 1]),
                                   x.count = x.count,
                                   rank.source = rank.source,
                                   comm = comm)
    } else{
      x.count <- spmd.bcast.integer(integer(COMM.SIZE),
                                    rank.source = rank.source, comm = comm)
      ret <- spmd.scatterv.integer(integer(0),
                                   integer(x.count[COMM.RANK + 1]),
                                   x.count = x.count,
                                   rank.source = rank.source,
                                   comm = comm)
    }
    dim(ret) <- ret.dim
    return(ret)
  }
} # End of spmd.scatter.array().

