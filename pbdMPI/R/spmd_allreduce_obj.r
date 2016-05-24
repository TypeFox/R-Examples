### Median level functions for R objects. These should not be in S3/S4.

### For general types.
spmd.allreduce.object <- function(x, op = .pbd_env$SPMD.CT$op,
    comm = .pbd_env$SPMD.CT$comm){
  x <- try(as.double(x), silent = TRUE)
  if(class(x) == "try-error"){
    stop(x)
  }
  .Call("spmd_allreduce_double", x, double(length(x)),
        which(op[1] == .pbd_env$SPMD.OP), as.integer(comm),
        PACKAGE = "pbdMPI")
} # End of spmd.allreduce.object().

### For array only.
spmd.allreduce.array <- function(x, op = .pbd_env$SPMD.CT$op,
    comm = .pbd_env$SPMD.CT$comm){
  COMM.SIZE <- spmd.comm.size(comm)

  all.check <- spmd.allreduce.integer(
                 as.integer(is.double(x) && length(x) > 0),
                 integer(1), op = "sum", comm = comm) == COMM.SIZE
  if(all.check){
    ret <- spmd.allreduce.double(x, double(length(x)), op = op[1], comm = comm)
    dim(ret) <- dim(x)
    return(ret)
  }

  all.check <- spmd.allreduce.integer(
                 as.integer(is.integer(x) && length(x) > 0),
                 integer(1), op = "sum", comm = comm) == COMM.SIZE
  if(all.check){
    ret <- spmd.allreduce.integer(x, integer(length(x)), op = op[1],
                                  comm = comm)
    dim(ret) <- dim(x)
    return(ret)
  }

  spmd.allreduce.object(x, op = op, comm = comm)
} # End of spmd.allreduce.array().

