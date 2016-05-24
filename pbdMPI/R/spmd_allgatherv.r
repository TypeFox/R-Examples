### S4 functions.

### Default method.
spmd.allgather.default <- function(x, x.buffer = NULL, x.count = NULL,
    displs = NULL, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  all.array <- spmd.allreduce.integer(as.integer(is.array(x)),
                                      integer(1), op = "sum",
                                      comm = comm) == spmd.comm.size(comm)
  if(all.array){
    spmd.allgather.array(x, comm = comm, unlist = unlist)
  } else{
    spmd.allgather.object(x, comm = comm, unlist = unlist)
  }
} # End of spmd.allgather.default().

spmd.allgatherv.default <- spmd.allgather.default


### For allgather and basic types.
spmd.allgather.integer <- function(x, x.buffer, x.count = NULL,
    displs = NULL, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  .Call("spmd_allgather_integer", x, x.buffer,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.allgather.double().

spmd.allgather.double <- function(x, x.buffer, x.count = NULL,
    displs = NULL, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  .Call("spmd_allgather_double", x, x.buffer,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.allgather.double().

spmd.allgather.raw <- function(x, x.buffer, x.count = NULL,
    displs = NULL, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  .Call("spmd_allgather_raw", x, x.buffer,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.allgather.raw().


### For allgatherv and basic types.
spmd.allgatherv.integer <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    comm = .pbd_env$SPMD.CT$comm, unlist = .pbd_env$SPMD.CT$unlist){
  .Call("spmd_allgatherv_integer", x, x.buffer, x.count, displs,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.allgatherv.integer().

spmd.allgatherv.double <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    comm = .pbd_env$SPMD.CT$comm, unlist = .pbd_env$SPMD.CT$unlist){
  .Call("spmd_allgatherv_double", x, x.buffer, x.count, displs,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.allgatherv.double().

spmd.allgatherv.raw <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    comm = .pbd_env$SPMD.CT$comm, unlist = .pbd_env$SPMD.CT$unlist){
  .Call("spmd_allgatherv_raw", x, x.buffer, x.count, displs,
        as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.allgatherv.raw().


### S4 methods.
setGeneric(
  name = "allgather",
  useAsDefault = spmd.allgather.default
)

### For allgather.
setMethod(
  f = "allgather",
  signature = signature(x = "ANY",
                        x.buffer = "missing",
                        x.count = "missing"),
  definition = spmd.allgather.default
)
setMethod(
  f = "allgather",
  signature = signature(x = "integer",
                        x.buffer = "integer",
                        x.count = "missing"),
  definition = spmd.allgather.integer
)
setMethod(
  f = "allgather",
  signature = signature(x = "numeric",
                        x.buffer = "numeric",
                        x.count = "missing"),
  definition = spmd.allgather.double
)
setMethod(
  f = "allgather",
  signature = signature(x = "raw",
                        x.buffer = "raw",
                        x.count = "missing"),
  definition = spmd.allgather.raw
)

### For allgatherv.
setMethod(
  f = "allgather",
  signature = signature(x = "ANY",
                        x.buffer = "missing",
                        x.count = "integer"),
  definition = spmd.allgatherv.default
)
setMethod(
  f = "allgather",
  signature = signature(x = "ANY",
                        x.buffer = "ANY",
                        x.count = "integer"),
  definition = spmd.allgatherv.default
)
setMethod(
  f = "allgather",
  signature = signature(x = "integer",
                        x.buffer = "integer",
                        x.count = "integer"),
  definition = spmd.allgatherv.integer
)
setMethod(
  f = "allgather",
  signature = signature(x = "numeric",
                        x.buffer = "numeric",
                        x.count = "integer"),
  definition = spmd.allgatherv.double
)
setMethod(
  f = "allgather",
  signature = signature(x = "raw",
                        x.buffer = "raw",
                        x.count = "integer"),
  definition = spmd.allgatherv.raw
)

