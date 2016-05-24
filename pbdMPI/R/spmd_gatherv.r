### S4 functions.

### Default method.
spmd.gather.default <- function(x, x.buffer = NULL, x.count = NULL,
    displs = NULL, rank.dest = .pbd_env$SPMD.CT$rank.root,
    comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  all.array <- spmd.allreduce.integer(as.integer(is.array(x)),
                                      integer(1), op = "sum",
                                      comm = comm) == spmd.comm.size(comm)
  if(all.array){
    spmd.gather.array(x, rank.dest = rank.dest, comm = comm, unlist = unlist)
  } else{
    spmd.gather.object(x, rank.dest = rank.dest, comm = comm, unlist = unlist)
  }
} # End of spmd.gather.default().

spmd.gatherv.default <- spmd.gather.default


### For gather and basic types.
spmd.gather.integer <- function(x, x.buffer, x.count = NULL, displs = NULL,
    rank.dest = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- .Call("spmd_gather_integer", x, x.buffer,
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.gather.double().

spmd.gather.double <- function(x, x.buffer, x.count = NULL, displs = NULL,
    rank.dest = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- .Call("spmd_gather_double", x, x.buffer,
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.gather.double().

spmd.gather.raw <- function(x, x.buffer, x.count = NULL, displs = NULL,
    rank.dest = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- .Call("spmd_gather_raw", x, x.buffer,
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.gather.raw().


### For gatherv and basic types.
spmd.gatherv.integer <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    rank.dest = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- .Call("spmd_gatherv_integer", x, x.buffer, x.count, displs,
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.gatherv.integer().

spmd.gatherv.double <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    rank.dest = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- .Call("spmd_gatherv_double", x, x.buffer, x.count, displs,
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.gatherv.double().

spmd.gatherv.raw <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    rank.dest = .pbd_env$SPMD.CT$rank.root, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- .Call("spmd_gatherv_raw", x, x.buffer, x.count, displs,
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.gatherv.raw().


### S4 methods.
setGeneric(
  name = "gather",
  useAsDefault = spmd.gather.default
)

### For gather.
setMethod(
  f = "gather",
  signature = signature(x = "ANY",
                        x.buffer = "missing",
                        x.count = "missing"),
  definition = spmd.gather.default
)
setMethod(
  f = "gather",
  signature = signature(x = "integer",
                        x.buffer = "integer",
                        x.count = "missing"),
  definition = spmd.gather.integer
)
setMethod(
  f = "gather",
  signature = signature(x = "numeric",
                        x.buffer = "numeric",
                        x.count = "missing"),
  definition = spmd.gather.double
)
setMethod(
  f = "gather",
  signature = signature(x = "raw",
                        x.buffer = "raw",
                        x.count = "missing"),
  definition = spmd.gather.raw
)

### For gatherv.
setMethod(
  f = "gather",
  signature = signature(x = "ANY",
                        x.buffer = "missing",
                        x.count = "integer"),
  definition = spmd.gatherv.default
)
setMethod(
  f = "gather",
  signature = signature(x = "ANY",
                        x.buffer = "ANY",
                        x.count = "integer"),
  definition = spmd.gatherv.default
)
setMethod(
  f = "gather",
  signature = signature(x = "integer",
                        x.buffer = "integer",
                        x.count = "integer"),
  definition = spmd.gatherv.integer
)
setMethod(
  f = "gather",
  signature = signature(x = "numeric",
                        x.buffer = "numeric",
                        x.count = "integer"),
  definition = spmd.gatherv.double
)
setMethod(
  f = "gather",
  signature = signature(x = "raw",
                        x.buffer = "raw",
                        x.count = "integer"),
  definition = spmd.gatherv.raw
)

