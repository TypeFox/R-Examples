### S4 functions.

### Default method.
spmd.irecv.default <- function(x.buffer = NULL,
    rank.source = .pbd_env$SPMD.CT$rank.source, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request,
    status = .pbd_env$SPMD.CT$status){
  ### This implementation is the same as spmd.recv.default(), because
  ### a blocking probe should be evoked to get the length of object for
  ### preparing enough x.buffer.
  spmd.recv.default(rank.source = rank.source, tag = tag, comm = comm,
                    status = status)
} # End of spmd.irecv.default().

### For irecv.
spmd.irecv.integer <- function(x.buffer,
    rank.source = .pbd_env$SPMD.CT$rank.source, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request,
    status = .pbd_env$SPMD.CT$status){
  .Call("spmd_irecv_integer", x.buffer, as.integer(rank.source),
        as.integer(tag), as.integer(comm), as.integer(request),
        PACKAGE = "pbdMPI")
} # End of spmd.irecv.integer().

spmd.irecv.double <- function(x.buffer,
    rank.source = .pbd_env$SPMD.CT$rank.source, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request,
    status = .pbd_env$SPMD.CT$status){
  .Call("spmd_irecv_double", x.buffer, as.integer(rank.source),
        as.integer(tag), as.integer(comm), as.integer(request),
        PACKAGE = "pbdMPI")
} # End of spmd.irecv.double().

spmd.irecv.raw <- function(x.buffer,
    rank.source = .pbd_env$SPMD.CT$rank.source, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request,
    status = .pbd_env$SPMD.CT$status){
  .Call("spmd_irecv_raw", x.buffer, as.integer(rank.source),
        as.integer(tag), as.integer(comm), as.integer(request),
        PACKAGE = "pbdMPI")
} # End of spmd.irecv.raw().


### S4 methods.
setGeneric(
  name = "irecv",
  useAsDefault = spmd.irecv.default
)

### For irecv.
setMethod(
  f = "irecv",
  signature = signature(x.buffer = "ANY"),
  definition = spmd.irecv.default
)
setMethod(
  f = "irecv",
  signature = signature(x.buffer = "integer"),
  definition = spmd.irecv.integer
)
setMethod(
  f = "irecv",
  signature = signature(x.buffer = "numeric"),
  definition = spmd.irecv.double
)
setMethod(
  f = "irecv",
  signature = signature(x.buffer = "raw"),
  definition = spmd.irecv.raw
)

