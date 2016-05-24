suppressMessages(library(pbdMPI, quietly = TRUE))

### Define method.
allgather.simple <- function(x, x.buffer = NULL, x.count = NULL,
    displs = NULL, comm = .pbd_env$SPMD.CT$comm,
    unlist = .pbd_env$SPMD.CT$unlist){
  ret <- allgather(as.integer(x@a[1]), x.buffer = integer(comm.size(comm)),
                   comm = comm)
  new("simple", a = ret)
}
.a <- setClass(Class = "simple",
               representation = representation(a = "numeric"))
.a <- setMethod(f = "allgather",
                signature = signature(x = "simple", x.buffer = "missing",
                                      x.count = "missing"),
                definition = allgather.simple)

### Run.
init()
x <- new("simple", a = comm.rank())
y <- allgather(x)
comm.print(y, all.rank = TRUE)
finalize()
