### S4 functions.

### Default method.
spmd.scatter.default <- function(x, x.buffer = NULL, x.count = NULL,
    displs = NULL, rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  if(spmd.comm.rank(comm) == rank.source){
    COMM.SIZE <- spmd.comm.size(comm)
    check <- c(is.list(x) && length(x) == spmd.comm.size(comm),
               do.call("sum", lapply(1:COMM.SIZE,
                                function(i) is.array(x[[i]]))) == COMM.SIZE)
    spmd.bcast.integer(as.integer(check), rank.source = rank.source,
                       comm = comm)
  } else{                             
    check <- spmd.bcast.integer(integer(2), rank.source = rank.source,
                                comm = comm)
  }

  if(!check[1]){
    stop("x should be a list and length COMM.SIZE")
  }
  if(check[2]){
    spmd.scatter.array(x, rank.source = rank.source, comm = comm)
  } else{              
    spmd.scatter.object(x, rank.source = rank.source, comm = comm)
  }
} # End of spmd.scatter.default().

spmd.scatterv.default <- spmd.scatter.default


### For scatter and basic types.
spmd.scatter.integer <- function(x, x.buffer, x.count = NULL, displs = NULL,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_scatter_integer", x, x.buffer,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.scatter.double().

spmd.scatter.double <- function(x, x.buffer, x.count = NULL, displs = NULL,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_scatter_double", x, x.buffer,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.scatter.double().

spmd.scatter.raw <- function(x, x.buffer, x.count = NULL, displs = NULL,
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_scatter_raw", x, x.buffer,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.scatter.raw().


### For scatterv and basic types.
spmd.scatterv.integer <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_scatterv_integer", x, x.buffer, x.count, displs,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.scatterv.integer().

spmd.scatterv.double <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_scatterv_double", x, x.buffer, x.count, displs,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.scatterv.double().

spmd.scatterv.raw <- function(x, x.buffer, x.count,
    displs = c(0L, cumsum(x.count)),
    rank.source = .pbd_env$SPMD.CT$rank.source,
    comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_scatterv_raw", x, x.buffer, x.count, displs,
        as.integer(rank.source), as.integer(comm), PACKAGE = "pbdMPI")
} # End of spmd.scatterv.raw().


### S4 methods.
setGeneric(
  name = "scatter",
  useAsDefault = spmd.scatter.default
)

### For scatter.
setMethod(
  f = "scatter",
  signature = signature(x = "ANY",
                        x.buffer = "missing",
                        x.count = "missing"),
  definition = spmd.scatter.default
)
setMethod(
  f = "scatter",
  signature = signature(x = "integer",
                        x.buffer = "integer",
                        x.count = "missing"),
  definition = spmd.scatter.integer
)
setMethod(
  f = "scatter",
  signature = signature(x = "numeric",
                        x.buffer = "numeric",
                        x.count = "missing"),
  definition = spmd.scatter.double
)
setMethod(
  f = "scatter",
  signature = signature(x = "raw",
                        x.buffer = "raw",
                        x.count = "missing"),
  definition = spmd.scatter.raw
)


### For scatterv.
setMethod(
  f = "scatter",
  signature = signature(x = "ANY",
                        x.buffer = "missing",
                        x.count = "integer"),
  definition = spmd.scatterv.default
)
setMethod(
  f = "scatter",
  signature = signature(x = "ANY",
                        x.buffer = "ANY",
                        x.count = "integer"),
  definition = spmd.scatterv.default
)
setMethod(
  f = "scatter",
  signature = signature(x = "integer",
                        x.buffer = "integer",
                        x.count = "integer"),
  definition = spmd.scatterv.integer
)
setMethod(
  f = "scatter",
  signature = signature(x = "numeric",
                        x.buffer = "numeric",
                        x.count = "integer"),
  definition = spmd.scatterv.double
)
setMethod(
  f = "scatter",
  signature = signature(x = "raw",
                        x.buffer = "raw",
                        x.count = "integer"),
  definition = spmd.scatterv.raw
)

