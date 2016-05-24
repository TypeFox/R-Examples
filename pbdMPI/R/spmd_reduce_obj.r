### Median level functions for R objects. These should not be in S3/S4.

### For general types.
spmd.reduce.object <- function(x, op = .pbd_env$SPMD.CT$op,
    rank.dest = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  x <- try(as.double(x), silent = TRUE)
  if(class(x) == "try-error"){
    stop(x) 
  }
  ret <- spmd.reduce.double(x, double(length(x)), op = op[1],
                            rank.dest = rank.dest, comm = comm)
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.reduce.object().

### For array only.
spmd.reduce.array <- function(x, op = .pbd_env$SPMD.CT$op,
    rank.dest = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- spmd.comm.size(comm)

  all.check <- spmd.allreduce.integer(
                 as.integer(is.double(x) && length(x) > 0),
                 integer(1), op = "sum", comm = comm) == COMM.SIZE
  if(all.check){
    ret <- spmd.reduce.double(x, double(length(x)),
                              op = op[1], rank.dest = rank.dest, comm = comm)
    if(spmd.comm.rank(comm) != rank.dest){
      return(invisible())
    }
    dim(ret) <- dim(x)
    return(ret)
  }

  all.check <- spmd.allreduce.integer(
                 as.integer(is.integer(x) && length(x) > 0),
                 integer(1), op = "sum", comm = comm) == COMM.SIZE
  if(all.check){
    ret <- spmd.reduce.integer(x, integer(length(x)),
                               op = op[1], rank.dest = rank.dest, comm = comm)
    if(spmd.comm.rank(comm) != rank.dest){
      return(invisible())
    }
    dim(ret) <- dim(x)
    return(ret)
  }

  spmd.reduce.object(x, op = op[1], rank.dest = rank.dest, comm = comm)
} # End of spmd.reduce.array().

