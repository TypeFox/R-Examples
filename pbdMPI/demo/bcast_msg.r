### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
x <- .comm.rank
msg <- NULL
if(.comm.rank == 0){
  msg <- "y <- as.character(x)"
}

### Get message.
msg <- spmd.bcast.message(msg)

### Evaluate the message: convert x to a string
eval(parse(text = msg))

z <- spmd.bcast.string(y)
comm.cat("\nBcast string:\n", quiet = TRUE)
comm.print(z, all.rank = TRUE)

### Finish.
finalize()
