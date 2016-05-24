### S4 functions.

### Default method.
spmd.reduce.default <- function(x, x.buffer = NULL, op = .pbd_env$SPMD.CT$op,
    rank.dest = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
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
    spmd.reduce.array(x, op = op[1], rank.dest = rank.dest, comm = comm)
  } else{
    spmd.reduce.object(x, op = op[1], rank.dest = rank.dest, comm = comm)
  }
} # End of spmd.reduce.default().


### For reduce and basic types.
spmd.reduce.integer <- function(x, x.buffer, op = .pbd_env$SPMD.CT$op,
    rank.dest = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  ret <- .Call("spmd_reduce_integer", x, x.buffer,
               which(op[1] == .pbd_env$SPMD.OP),
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.reduce.integer().

spmd.reduce.double <- function(x, x.buffer, op = .pbd_env$SPMD.CT$op,
    rank.dest = .pbd_env$SPMD.CT$rank.source, comm = .pbd_env$SPMD.CT$comm){
  ret <- .Call("spmd_reduce_double", x, x.buffer,
               which(op[1] == .pbd_env$SPMD.OP),
               as.integer(rank.dest), as.integer(comm), PACKAGE = "pbdMPI")
  if(spmd.comm.rank(comm) != rank.dest){
    return(invisible())
  }
  ret
} # End of spmd.reduce.double().


### S4 methods.
setGeneric(
  name = "reduce",
  useAsDefault = spmd.reduce.default
)

### For reduce.
setMethod(
  f = "reduce",
  signature = signature(x = "ANY", x.buffer = "missing"),
  definition = spmd.reduce.default
)
setMethod(
  f = "reduce",
  signature = signature(x = "integer", x.buffer = "integer"),
  definition = spmd.reduce.integer
)
setMethod(
  f = "reduce",
  signature = signature(x = "numeric", x.buffer = "numeric"),
  definition = spmd.reduce.double
)
