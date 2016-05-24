### S4 functions.

### Default method.
spmd.bcast.default <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) == rank.source){
    is.check <- spmd.bcast.integer(as.integer(is.array(x)), comm = comm)
  } else{
    is.check <- spmd.bcast.integer(integer(1), comm = comm)
  }

  if(is.check[1]){
    ret <- spmd.bcast.array(x, rank.source = rank.source, comm = comm)
  } else{
    ret <- spmd.bcast.object(x, rank.source = rank.source, comm = comm)
  }

  ret
} # End of spmd.bcast.default().

### For bcast and basic types.
spmd.bcast.integer <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_bcast_integer", x,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.bcast.integer().

spmd.bcast.double <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_bcast_double", x,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.bcast.double().

spmd.bcast.raw <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_bcast_raw", x,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.bcast.raw().

spmd.bcast.string <- function(x,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_bcast_string", x[1],
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.bcast.string().


### S4 methods.
setGeneric(
  name = "bcast",
  useAsDefault = spmd.bcast.default
)

### For bcast.
setMethod(
  f = "bcast",
  signature = signature(x = "ANY"),
  definition = spmd.bcast.default
)
setMethod(
  f = "bcast",
  signature = signature(x = "integer"),
  definition = spmd.bcast.default
)
setMethod(
  f = "bcast",
  signature = signature(x = "numeric"),
  definition = spmd.bcast.default
)
setMethod(
  f = "bcast",
  signature = signature(x = "raw"),
  definition = spmd.bcast.raw
)

