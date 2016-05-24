### S3 tool function.

comm.any <- function(x, na.rm = FALSE, comm = .pbd_env$SPMD.CT$comm){
  ret <- spmd.allgather.integer(any(x, na.rm = na.rm),
                                integer(comm.size(comm)), comm = comm)
  any(ret, na.rm = na.rm)
} # End of comm.any().

comm.all <- function(x, na.rm = FALSE, comm = .pbd_env$SPMD.CT$comm){
  ret <- spmd.allgather.integer(all(x, na.rm = na.rm),
                                integer(comm.size(comm)), comm = comm)
  all(ret, na.rm = na.rm)
} # End of comm.all().
