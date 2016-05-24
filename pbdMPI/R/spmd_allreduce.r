### S4 functions.

### Default method.
spmd.allreduce.default <- function(x, x.buffer = NULL,
    op = .pbd_env$SPMD.CT$op, comm = .pbd_env$SPMD.CT$comm){
  op <- match.arg(tolower(op[1]), .pbd_env$SPMD.OP)

  if(op %in% c("land", "band", "lor", "bor", "lxor", "bxor")){
    x <- as.integer(x)
    if(!is.null(x.buffer)){
      x.buffer <- as.integer(x.buffer)
    }
  }

  all.array <- spmd.allreduce.integer(
                   as.integer(is.array(x) && length(x) > 0),
                   integer(1), op = "sum",
                   comm = comm) == spmd.comm.size(comm)
  if(all.array){
    x <- spmd.allcheck.type(x, comm = comm)
    spmd.allreduce.array(x, op = op[1], comm = comm)
  } else{
    spmd.allreduce.object(x, op = op[1], comm = comm)
  }
} # End of spmd.allreduce.default().


### For allreduce and basic types.
spmd.allreduce.integer <- function(x, x.buffer,
    op = .pbd_env$SPMD.CT$op, comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_allreduce_integer", x, x.buffer,
        which(op[1] == .pbd_env$SPMD.OP), as.integer(comm),
        PACKAGE = "pbdMPI")
} # End of spmd.allreduce.integer().

spmd.allreduce.double <- function(x, x.buffer,
    op = .pbd_env$SPMD.CT$op, comm = .pbd_env$SPMD.CT$comm){
  .Call("spmd_allreduce_double", x, x.buffer,
        which(op[1] == .pbd_env$SPMD.OP), as.integer(comm),
        PACKAGE = "pbdMPI")
} # End of spmd.allreduce.double().


### S4 methods.
setGeneric(
  name = "allreduce",
  useAsDefault = spmd.allreduce.default
)

### For allreduce.
setMethod(
  f = "allreduce",
  signature = signature(x = "ANY", x.buffer = "missing"),
  definition = spmd.allreduce.default
)
setMethod(
  f = "allreduce",
  signature = signature(x = "integer", x.buffer = "integer"),
  definition = spmd.allreduce.integer
)
setMethod(
  f = "allreduce",
  signature = signature(x = "numeric", x.buffer = "numeric"),
  definition = spmd.allreduce.double
)

