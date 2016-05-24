### Utilities.

spmd.probe <- function(rank.source = .pbd_env$SPMD.CT$rank.source,
    tag = .pbd_env$SPMD.CT$tag, comm = .pbd_env$SPMD.CT$comm,
    status = .pbd_env$SPMD.CT$status){
  ret <- .Call("spmd_probe", as.integer(rank.source), as.integer(tag),
               as.integer(comm), as.integer(status), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.probe().

probe <- spmd.probe

spmd.iprobe <- function(rank.source = .pbd_env$SPMD.CT$rank.source,
    tag = .pbd_env$SPMD.CT$tag, comm = .pbd_env$SPMD.CT$comm,
    status = .pbd_env$SPMD.CT$status){
  ret <- .Call("spmd_iprobe", as.integer(rank.source), as.integer(tag),
               as.integer(comm), as.integer(status), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.iprobe().

iprobe <- spmd.iprobe

spmd.anysource <- function(){
  .Call("spmd_anysource", PACKAGE = "pbdMPI")
} # End of spmd.anysource().

anysource <- spmd.anysource

spmd.anytag <- function(){
  .Call("spmd_anytag", PACKAGE = "pbdMPI")
} # End of spmd.anytag().

anytag <- spmd.anytag

spmd.get.sourcetag <- function(status = .pbd_env$SPMD.CT$status){
  .Call("spmd_get_sourcetag", as.integer(status), PACKAGE = "pbdMPI")
} # End of spmd.get.sourcetag().

get.sourcetag <- spmd.get.sourcetag

spmd.get.count <- function(data.type, status = .pbd_env$SPMD.CT$status){
  .Call("spmd_get_count", as.integer(data.type), as.integer(status),
        PACKAGE = "pbdMPI")
} # End of spmd.get.count().

get.count <- spmd.get.count

spmd.is.comm.null <- function(comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_is_comm_null", as.integer(comm),
        PACKAGE = "pbdMPI")
} # End of spmd.is.comm.null().

is.comm.null <- spmd.is.comm.null
