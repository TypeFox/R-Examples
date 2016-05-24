### S4 functions.

### Default method.
spmd.isend.default <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request){
  ### WCC: This isend() should go with wait(), otherwise the new R object,
  ###      "serialize(x, NULL)", is NOT sent correctly since it is not
  ###      protected by R when the call "spmd.isend.default()" is returned.
  ###      Add wait() is equivalent to use send() function.
  # spmd.isend.raw(serialize(x, NULL), rank.dest = as.integer(rank.dest),
  #                tag = as.integer(tag), comm = as.integer(comm),
  #                request = as.integer(request))
  ### This implementation is the same as spmd.send.default(), because
  ### a blocking wait should be evoked to make sure buffer is sent completely.
  ### But this can cause dead lock.
  # spmd.send.raw(serialize(x, NULL), rank.dest = rank.dest,
  #               tag = tag, comm = comm)
  ### Use non-blocking buffer to avoid dead lock and use non-block send.
  if(is.null(.pbd_env$SPMD.NB.BUFFER)){
    .pbd_env$SPMD.NB.BUFFER <- list()
  }
  .pbd_env$SPMD.NB.BUFFER[[length(.pbd_env$SPMD.NB.BUFFER) + 1]] <- serialize(x, NULL)
  spmd.isend.raw(.pbd_env$SPMD.NB.BUFFER[[length(.pbd_env$SPMD.NB.BUFFER)]],
                 rank.dest = as.integer(rank.dest),
                 tag = as.integer(tag), comm = as.integer(comm),
                 request = as.integer(request))
  invisible()
} # End of spmd.isend.default().

### For isend.
spmd.isend.integer <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request){
  .Call("spmd_isend_integer", x, as.integer(rank.dest), as.integer(tag),
        as.integer(comm), as.integer(request), PACKAGE = "pbdMPI")
  invisible()
} # End of spmd.isend.integer().

spmd.isend.double <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request){
  .Call("spmd_isend_double", x, as.integer(rank.dest), as.integer(tag),
        as.integer(comm), as.integer(request), PACKAGE = "pbdMPI")
  invisible()
} # End of spmd.isend.double().

spmd.isend.raw <- function(x,
    rank.dest = .pbd_env$SPMD.CT$rank.dest, tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, request = .pbd_env$SPMD.CT$request){
  .Call("spmd_isend_raw", x, as.integer(rank.dest), as.integer(tag),
        as.integer(comm), as.integer(request), PACKAGE = "pbdMPI")
  invisible()
} # End of spmd.isend.double().


### S4 methods.
setGeneric(
  name = "isend",
  useAsDefault = spmd.isend.default
)

### For isend.
setMethod(
  f = "isend",
  signature = signature(x = "ANY"),
  definition = spmd.isend.default
)
setMethod(
  f = "isend",
  signature = signature(x = "integer"),
  definition = spmd.isend.integer
)
setMethod(
  f = "isend",
  signature = signature(x = "numeric"),
  definition = spmd.isend.double
)
setMethod(
  f = "isend",
  signature = signature(x = "raw"),
  definition = spmd.isend.raw
)

