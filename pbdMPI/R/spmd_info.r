### Info functions.

spmd.info.create <- function(info = .pbd_env$SPMD.CT$info){
  ret <- .Call("spmd_info_create", as.integer(info), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.info.create().

info.create <- spmd.info.create

spmd.info.set <- function(info = .pbd_env$SPMD.CT$info, key, value){
  ret <- .Call("spmd_info_set", as.integer(info), as.character(key),
               as.character(value), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.info.set().

info.set <- spmd.info.set

spmd.info.free <- function(info = .pbd_env$SPMD.CT$info){
  ret <- .Call("spmd_info_free", as.integer(info), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.info.free().

info.free <- spmd.info.free


### Fortran supporting function.
spmd.info.c2f <- function(info = .pbd_env$SPMD.CT$info){
  .Call("spmd_info_c2f", as.integer(info), PACKAGE = "pbdMPI")
} # End of spmd.info.c2f().

info.c2f <- spmd.info.c2f
