### All common functions.

comm.allcommon <- function(x, comm = .pbd_env$SPMD.CT$comm,
    lazy.check = .pbd_env$SPMD.CT$lazy.check){
  if(lazy.check){
    ### Faster but dangerous.
    ret <- allreduce(x, comm = comm) / spmd.comm.size(comm) == x
  } else{ 
    ### Much slower but safer.
    tmp <- do.call("cbind", allgather(x, comm = comm))
    ret <- apply(tmp, 1, function(x){ length(unique(x)) }) == 1
  }
  ret 
} # End of comm.allcommon().


### Internal only.
comm.allcommon.integer <- function(x, comm = .pbd_env$SPMD.CT$comm){
  x <- as.integer(x)
  spmd.allreduce.integer(x, integer(length(x)), comm = comm, op = "band") == x
} # End of comm.allcommon.integer().

