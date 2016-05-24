### S4 functions.

### Default method.
spmd.sendrecv.replace.default <- function(x,
    rank.dest = (comm.rank(.pbd_env$SPMD.CT$comm) + 1) %%
                comm.size(.pbd_env$SPMD.CT$comm),
    send.tag = .pbd_env$SPMD.CT$tag,
    rank.source = (comm.rank(.pbd_env$SPMD.CT$comm) - 1) %%
                  comm.size(.pbd_env$SPMD.CT$comm),
    recv.tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  x.raw <- serialize(x, NULL)
  total.org <- length(x.raw)
  total.new <- spmd.sendrecv.replace.integer(as.integer(total.org),
                 rank.dest = rank.source,
                 send.tag = send.tag,
                 rank.source = rank.dest,
                 recv.tag = recv.tag,
                 comm = comm, status = status)
  if(total.org == total.new){
    unserialize(spmd.sendrecv.replace.raw(x.raw,
                  rank.dest = rank.dest,
                  send.tag = send.tag,
                  rank.source = rank.source,
                  recv.tag = recv.tag,
                  comm = comm, status = status))
  } else{
    stop("Objects are not consistent.")
  }
} # End of spmd.sendrecv.replace.default().


### For sendrecv.replace.
spmd.sendrecv.replace.integer <- function(x,
    rank.dest = (comm.rank(.pbd_env$SPMD.CT$comm) + 1) %%
                comm.size(.pbd_env$SPMD.CT$comm),
    send.tag = .pbd_env$SPMD.CT$tag,
    rank.source = (comm.rank(.pbd_env$SPMD.CT$comm) - 1) %%
                  comm.size(.pbd_env$SPMD.CT$comm),
    recv.tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  .Call("spmd_sendrecv_replace_integer", x,
        as.integer(rank.dest), as.integer(send.tag), as.integer(rank.source),
        as.integer(recv.tag),
        as.integer(comm), as.integer(status), PACKAGE = "pbdMPI")
} # End of spmd.sendrecv.replace.integer().

spmd.sendrecv.replace.double <- function(x,
    rank.dest = (comm.rank(.pbd_env$SPMD.CT$comm) + 1) %%
                comm.size(.pbd_env$SPMD.CT$comm),
    send.tag = .pbd_env$SPMD.CT$tag,
    rank.source = (comm.rank(.pbd_env$SPMD.CT$comm) - 1) %%
                  comm.size(.pbd_env$SPMD.CT$comm),
    recv.tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  .Call("spmd_sendrecv_replace_double", x,
        as.integer(rank.dest), as.integer(send.tag), as.integer(rank.source),
        as.integer(recv.tag),
        as.integer(comm), as.integer(status), PACKAGE = "pbdMPI")
} # End of spmd.sendrecv.replace.double().

spmd.sendrecv.replace.raw <- function(x,
    rank.dest = (comm.rank(.pbd_env$SPMD.CT$comm) + 1) %%
                comm.size(.pbd_env$SPMD.CT$comm),
    send.tag = .pbd_env$SPMD.CT$tag,
    rank.source = (comm.rank(.pbd_env$SPMD.CT$comm) - 1) %%
                  comm.size(.pbd_env$SPMD.CT$comm),
    recv.tag = .pbd_env$SPMD.CT$tag,
    comm = .pbd_env$SPMD.CT$comm, status = .pbd_env$SPMD.CT$status){
  .Call("spmd_sendrecv_replace_raw", x,
        as.integer(rank.dest), as.integer(send.tag), as.integer(rank.source),
        as.integer(recv.tag),
        as.integer(comm), as.integer(status), PACKAGE = "pbdMPI")
} # End of spmd.sendrecv.replace.raw().


### S4 methods for recv.
setGeneric(
  name = "sendrecv.replace",
  useAsDefault = spmd.sendrecv.replace.default
)
setMethod(
  f = "sendrecv.replace",
  signature = signature(x = "ANY"),
  definition = spmd.sendrecv.replace.default
)
setMethod(
  f = "sendrecv.replace",
  signature = signature(x = "integer"),
  definition = spmd.sendrecv.replace.integer
)
setMethod(
  f = "sendrecv.replace",
  signature = signature(x = "numeric"),
  definition = spmd.sendrecv.replace.double
)
setMethod(
  f = "sendrecv.replace",
  signature = signature(x = "raw"),
  definition = spmd.sendrecv.replace.raw
)

