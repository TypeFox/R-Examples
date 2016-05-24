# For computing global range.

comm.range <- function(..., na.rm = FALSE, comm = .pbd_env$SPMD.CT$comm){
  ret <- range(...)
  ret[1] <- min(do.call("c", allgather(ret[1], comm = comm)), na.rm = na.rm)
  ret[2] <- max(do.call("c", allgather(ret[2], comm = comm)), na.rm = na.rm)
  ret
} # End of comm.range().

comm.max <- function(..., na.rm = FALSE, comm = .pbd_env$SPMD.CT$comm){
  ret <- max(..., na.rm = na.rm)
  ret <- max(do.call("c", allgather(ret[1], comm = comm)), na.rm = na.rm)
  ret
} # End of comm.max().

comm.min <- function(..., na.rm = FALSE, comm = .pbd_env$SPMD.CT$comm){
  ret <- min(..., na.rm = na.rm)
  ret <- min(do.call("c", allgather(ret[1], comm = comm)), na.rm = na.rm)
  ret
} # End of comm.min().

