spmd.barrier <- function(comm = .pbd_env$SPMD.CT$comm){
  ret <- .Call("spmd_barrier", as.integer(comm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.barrier().

barrier <- spmd.barrier

spmd.comm.set.errhandler <- function(comm = .pbd_env$SPMD.CT$comm){
  ret <- .Call("spmd_comm_set_errhandler", as.integer(comm),
               PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.set.errhandler().

spmd.comm.is.null <- function(comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_comm_is_null", as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.comm.is.null().

comm.is.null <- spmd.comm.is.null

spmd.comm.rank <- function(comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_comm_rank", as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.comm.rank().

comm.rank <- spmd.comm.rank

spmd.comm.size <- function(comm = .pbd_env$SPMD.CT$comm){
  tmp <- .Call("spmd_comm_is_null", as.integer(comm), PACKAGE = "pbdMPI")

  if(tmp == 1){
    0L
  } else{
    .Call("spmd_comm_size", as.integer(comm), PACKAGE = "pbdMPI")
  }
} # End of spmd.comm.size().

comm.size <- spmd.comm.size

spmd.comm.dup <- function(comm, newcomm){
  ret <- .Call("spmd_comm_dup", as.integer(comm), as.integer(newcomm),
               PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.dup().

comm.dup <- spmd.comm.dup

spmd.comm.free <- function(comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.size(comm) == 0){
    stop(paste("It seems no members (workers) associated with comm", comm))
  }
  ret <- .Call("spmd_comm_free", as.integer(comm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.free().

comm.free <- spmd.comm.free

spmd.init <- function(set.seed = TRUE){
  ### Check even ".__DISABLE_MPI_INIT" is set by external API.
  # if(! exists(".__DISABLE_MPI_INIT__", envir = .GlobalEnv) ||
  #    get(".__DISABLE_MPI_INIT__", envir = .GlobalEnv) != TRUE){
  #   assign(".__DISABLE_MPI_INIT__", FALSE, envir = .GlobalEnv)
  # }

  ### We still need to initial memory for our own communicators.
  ### Copy the COMM_WORLD to the comm 0.
  ret <- .Call("spmd_initialize", PACKAGE = "pbdMPI")
  # assign(".comm.size", spmd.comm.size(), envir = .GlobalEnv)
  # assign(".comm.rank", spmd.comm.rank(), envir = .GlobalEnv)

  ### For seed.
  if(set.seed){
   # seed <- as.integer(Sys.getpid() + Sys.time())
    seed <- as.integer(runif(6, 1L, 2147483647L))
    seed <- .Call("spmd_bcast_integer", seed, 0L, 0L, PACKAGE = "pbdMPI")
    # seed <- rep(seed, 6)

    comm.size <- .Call("spmd_comm_size", 0L, PACKAGE = "pbdMPI")
    comm.rank <- .Call("spmd_comm_rank", 0L, PACKAGE = "pbdMPI")
    names <- as.character(0:(comm.size - 1))
    name <- as.character(comm.rank)

    suppressWarnings(eval(.lec.old.kind <- RNGkind(), envir = .GlobalEnv))
    suppressWarnings(eval(.lec.SetPackageSeed(seed), envir = .GlobalEnv))
    suppressWarnings(eval(.lec.CreateStream(names), envir = .GlobalEnv))
    suppressWarnings(eval(.lec.CurrentStream(name), envir = .GlobalEnv))
  }

  invisible(ret)
} # End of spmd.init().

init <- spmd.init

spmd.finalize <- function(mpi.finalize = .pbd_env$SPMD.CT$mpi.finalize){
  ### Do not remove ".__DISABLE_MPI_INIT__", leave it in .GlobalEnv for later
  ### uses.

  ### Only free the memory. Manually shut down MPI by "mpi.finalize".
  ### Let users take care of MPI shut down business.
  ret <- .Call("spmd_finalize", mpi.finalize, PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.finalize().

finalize <- spmd.finalize

spmd.is.finalized <- function(){
  ret <- .Call("spmd_is_finalized", PACKAGE = "pbdMPI")
  invisible(as.logical(ret))
} # End of spmd.finalized().

is.finalized <- spmd.is.finalized

spmd.is.master <- function(){
  tmp <- is.loaded("spmd_comm_get_parent", PACKAGE = "pbdMPI")
  if(tmp){
    as.logical(.Call("spmd_is_master", PACKAGE = "pbdMPI"))
  } else{
    if(spmd.comm.size(1L) > 0){
      spmd.comm.rank(1L) == 0
    } else{
      spmd.comm.rank(0L) == 0
    }
  }
} # End of spmd.is.master().

is.master <- spmd.is.master

spmd.get.processor.name <- function(short = TRUE){
  name <- .Call("spmd_get_processor_name", PACKAGE = "pbdMPI")
  if(short){
    name <- unlist(strsplit(name, "\\."))[1]
  }
  name
} # End of spmd.get.processor.name().

get.processor.name <- spmd.get.processor.name

spmd.comm.abort <- function(errorcode = 1, comm = .pbd_env$SPMD.CT$comm){
  ret <- .Call("spmd_comm_abort", as.integer(comm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.abort().

comm.abort <- spmd.comm.abort

spmd.comm.split <- function(comm = .pbd_env$SPMD.CT$comm, color = 0L,
    key = 0L, newcomm = .pbd_env$SPMD.CT$newcomm){
  ret <- .Call("spmd_comm_split", as.integer(comm), as.integer(color),
               as.integer(key), as.integer(newcomm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.split().

comm.split <- spmd.comm.split

spmd.comm.disconnect <- function(comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.size(comm)== 0){
    stop(paste("It seems no members (workers) associated with comm", comm))
  }
  ret <- .Call("spmd_comm_disconnect", as.integer(comm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.disconnect().

comm.disconnect <- spmd.comm.disconnect

spmd.comm.connect <- function(port.name,
    info = .pbd_env$SPMD.CT$info, rank.root = .pbd_env$SPMD.CT$rank.root,
    comm = .pbd_env$SPMD.CT$comm, newcomm = .pbd_env$SPMD.CT$newcomm){
  ret <- .Call("spmd_comm_connect", as.character(port.name),
               as.integer(info), as.integer(rank.root),
               as.integer(comm), as.integer(newcomm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.comm.connect().

comm.connect <- spmd.comm.connect

spmd.comm.accept <- function(port.name,
    info = .pbd_env$SPMD.CT$info, rank.root = .pbd_env$SPMD.CT$rank.root,
    comm = .pbd_env$SPMD.CT$comm, newcomm = .pbd_env$SPMD.CT$newcomm){
  ret <- .Call("spmd_comm_accept", as.character(port.name),
               as.integer(info), as.integer(rank.root),
               as.integer(comm), as.integer(newcomm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End spmd.comm.accept().

comm.accept <- spmd.comm.accept

spmd.port.open <- function(info = .pbd_env$SPMD.CT$info){
  port.name <- .Call("spmd_port_open", as.integer(info), PACKAGE = "pbdMPI")
  port.name
} # End spmd.port.open().

port.open <- spmd.port.open

spmd.port.close <- function(port.name){
  ret <- .Call("spmd_port_close", as.character(port.name), PACKAGE = "pbdMPI")
  invisible(ret)
} # End spmd.port.close().

port.close <- spmd.port.close

spmd.serv.publish <- function(port.name,
    serv.name = .pbd_env$SPMD.CT$serv.name,
    info = .pbd_env$SPMD.CT$info){
  ret <- .Call("spmd_serv_publish", as.character(serv.name),
               as.integer(info), as.character(port.name), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.serv.publish().

serv.publish <- spmd.serv.publish

spmd.serv.unpublish <- function(port.name,
    serv.name = .pbd_env$SPMD.CT$serv.name,
    info = .pbd_env$SPMD.CT$info){
  ret <- .Call("spmd_serv_unpublish", as.character(serv.name),
               as.integer(info), as.character(port.name), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.serv.unpublish().

serv.unpublish <- spmd.serv.unpublish

spmd.serv.lookup <- function(serv.name = .pbd_env$SPMD.CT$serv.name,
    info = .pbd_env$SPMD.CT$info){
  port.name <- .Call("spmd_serv_lookup", as.character(serv.name),
                     as.integer(info), PACKAGE = "pbdMPI")
  port.name
} # End of spmd.serv.lookup().

serv.lookup <- spmd.serv.lookup

spmd.comm.get.parent <- function(comm = .pbd_env$SPMD.CT$intercomm){
  .Call("spmd_comm_get_parent", as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.comm.get.parent().

spmd.intercomm.merge <- function(intercomm = .pbd_env$SPMD.CT$intercomm,
    high = 0L, comm = .pbd_env$SPMD.CT$comm){
  ret <- .Call("spmd_intercomm_merge", as.integer(intercomm), as.integer(high),
               as.integer(comm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.intercomm.merge().

intercomm.merge <- spmd.intercomm.merge

spmd.intercomm.create <- function(local.comm = .pbd_env$SPMD.CT$comm,
    local.leader = .pbd_env$SPMD.CT$rank.source,
    peer.comm = .pbd_env$SPMD.CT$intercomm,
    remote.leader = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    newintercomm = .pbd_env$SPMD.CT$newcomm){
  ret <- .Call("spmd_intercomm_create", as.integer(local.comm),
               as.integer(local.leader), as.integer(peer.comm),
               as.integer(remote.leader), as.integer(tag),
               as.integer(newintercomm), PACKAGE = "pbdMPI")
  invisible(ret)
} # End of spmd.intercomm.merge().

intercomm.create <- spmd.intercomm.create


### Fortran supporting function.
spmd.comm.c2f <- function(comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_comm_c2f", as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.comm.c2f().

comm.c2f <- spmd.comm.c2f
