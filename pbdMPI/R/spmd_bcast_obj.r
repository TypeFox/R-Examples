### Median level functions for R objects. These should not be in S3/S4.

### For general types.
spmd.bcast.object <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) == rank.source){
    x.raw <- serialize(x, NULL)
    spmd.bcast.integer(length(x.raw), rank.source = rank.source, comm = comm)
    spmd.bcast.raw(x.raw, rank.source = rank.source, comm = comm)
    return(x) 
  } else{
    x.count <- spmd.bcast.integer(integer(1), rank.source = rank.source,
                                  comm = comm)
    unserialize(spmd.bcast.raw(raw(x.count), rank.source = rank.source,
                               comm = comm))
  }
} # End of spmd.bcast.object().

### For array only.
spmd.bcast.array <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) == rank.source){
    spmd.bcast.integer(length(dim(x)), rank.source = rank.source, comm = comm)
    spmd.bcast.integer(dim(x), rank.source = rank.source, comm = comm)

    check <- spmd.bcast.integer(as.integer(is.double(x)),
                                rank.source = rank.source, comm = comm)
    if(check){
      spmd.bcast.double(x, rank.source = rank.source, comm = comm)
      return(x)
    }

    check <- spmd.bcast.integer(as.integer(is.integer(x)),
                                rank.source = rank.source, comm = comm)
    if(check){
      spmd.bcast.integer(x, rank.source = rank.source, comm = comm)
      return(x)
    }
  } else{
    n.dim <- spmd.bcast.integer(integer(1),
                                rank.source = rank.source, comm = comm)
    x.dim <- spmd.bcast.integer(integer(n.dim),
                                rank.source = rank.source, comm = comm)

    check <- spmd.bcast.integer(integer(1),
                                rank.source = rank.source, comm = comm)
    if(check){
      ret <- spmd.bcast.double(double(prod(x.dim)),
                               rank.source = rank.source, comm = comm)
      dim(ret) <- x.dim
      return(ret)
    }

    check <- spmd.bcast.integer(integer(1),
                                rank.source = rank.source, comm = comm)
    if(check){
      ret <- spmd.bcast.integer(integer(prod(x.dim)),
                               rank.source = rank.source, comm = comm)
      dim(ret) <- x.dim
      return(ret)
    }
  }

  spmd.bcast.object(x, rank.source = rank.source, comm = comm) 
} # End of spmd.bcast.array().

### For message.
spmd.bcast.message <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) == rank.source){
    spmd.bcast.integer(nchar(x[1]), rank.source = rank.source, comm = comm)
    spmd.bcast.string(x[1], rank.source = rank.source, comm = comm)
    return(x[1])
  } else{
    x.count <- spmd.bcast.integer(integer(1), rank.source = rank.source,
                                  comm = comm)
    spmd.bcast.string(paste0(rep(" ", x.count), collapse = ""),
                      rank.source = rank.source, comm = comm)
  }
} # End of spmd.bcast.message().
