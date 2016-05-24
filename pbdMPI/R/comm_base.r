### Some base functions.

comm.length <- function(x, comm = .pbd_env$SPMD.CT$comm){
  if(! comm.all(is.vector(x), comm = comm)){
    comm.stop("x should be a vector.", comm = comm)
  }

  ret <- allreduce(length(x), op = "sum", comm = comm)
  ret
} # End of comm.length().

comm.sum <- function(..., na.rm = TRUE, comm = .pbd_env$SPMD.CT$comm){
  tmp <- sum(..., na.rm = na.rm)
  ret <- allreduce(tmp, op = "sum", comm = comm)
  ret
} # End of comm.sum().

comm.mean <- function(x, na.rm = TRUE, comm = .pbd_env$SPMD.CT$comm){
  if(! comm.all(is.vector(x), comm = comm)){
    comm.stop("x should be a vector.", comm = comm)
  }

  if(na.rm){
    id.ok <- (!is.na(x)) & (!is.nan(x))
    x <- x[id.ok]
  }
  tl.x <- allreduce(length(x), op = "sum", comm = comm)

  ret <- comm.sum(x / tl.x, comm = comm)
  ret
} # End of comm.mean().

comm.var <- function(x, na.rm = TRUE, comm = .pbd_env$SPMD.CT$comm){
  if(! comm.all(is.vector(x), comm = comm)){
    comm.stop("x should be a vector.", comm = comm)
  }

  if(na.rm){
    id.ok <- (!is.na(x)) & (!is.nan(x))
    x <- x[id.ok]
  }
  tl.x <- allreduce(length(x), op = "sum", comm = comm)

  mean.x.2 <- comm.sum(x^2 / (tl.x - 1), comm = comm)
  mean.x <- comm.sum(x / tl.x, comm = comm)
  ret <- mean.x.2 - mean.x^2 * (tl.x / (tl.x - 1))
  ret
} # End of comm.var().

comm.sd <- function(x, na.rm = TRUE, comm = .pbd_env$SPMD.CT$comm){
  ret <- comm.var(x, na.rm = na.rm, comm = comm)
  ret <- sqrt(ret)
  ret
} # End of comm.sd().

