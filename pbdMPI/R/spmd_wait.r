### For non-blocking calls.

spmd.wait <- function(request = .pbd_env$SPMD.CT$request,
    status = .pbd_env$SPMD.CT$status){
  ret <- .Call("spmd_wait", as.integer(request), as.integer(status),
               PACKAGE = "pbdMPI")
  ### Clear non-blocking buffer.
  .pbd_env$SPMD.NB.BUFFER <- list()
  invisible(ret)
} # End of spmd.wait().

wait <- spmd.wait

spmd.waitany <- function(count, status = .pbd_env$SPMD.CT$status){
  ret <- .Call("spmd_waitany",  as.integer(count), as.integer(status),
               PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.waitany().

waitany <- spmd.waitany

spmd.waitsome <- function(count){
  tmp <- .Call("spmd_waitsome",  as.integer(count), PACKAGE = "pbdMPI")
  if(tmp[1] < 0 || tmp[1] > count){
    return(list(count = tmp[1], indices = NULL))
  } else{
    return(list(count = tmp[1], indices = tmp[2:(1 + tmp[1])]))
  }
} # End of spmd.waitsome().

waitsome <- spmd.waitsome

spmd.waitall <- function(count){
  ret <- .Call("spmd_waitall", as.integer(count), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.waitall().

waitall <- spmd.waitall
