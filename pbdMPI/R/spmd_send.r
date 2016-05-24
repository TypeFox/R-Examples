### S4 functions.

### Default method.
spmd.send.default <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm){
  spmd.send.raw(serialize(x, NULL), rank.dest = rank.dest,
                tag = tag, comm = comm)
  invisible()
} # End of spmd.send.default().


### For send.
spmd.send.integer <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_send_integer", x, as.integer(rank.dest), as.integer(tag),
        as.integer(comm), PACKAGE = "pbdMPI")
  invisible()
} # End of spmd.send.integer().

spmd.send.double <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_send_double", x, as.integer(rank.dest), as.integer(tag),
        as.integer(comm), PACKAGE = "pbdMPI")
  invisible()
} # End of spmd.send.double().

spmd.send.raw <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_send_raw", x, as.integer(rank.dest), as.integer(tag),
        as.integer(comm), PACKAGE = "pbdMPI")
  invisible()
} # End of spmd.send.raw().


### S4 methods.
setGeneric(
  name = "send",
  useAsDefault = spmd.send.default
)

### For send.
setMethod(
  f = "send",
  signature = signature(x = "ANY"),
  definition = spmd.send.default
)
setMethod(
  f = "send",
  signature = signature(x = "integer"),
  definition = spmd.send.integer
)
setMethod(
  f = "send",
  signature = signature(x = "numeric"),
  definition = spmd.send.double
)
setMethod(
  f = "send",
  signature = signature(x = "raw"),
  definition = spmd.send.raw
)

