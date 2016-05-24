### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
n <- as.integer(2)
x <- 1:(.comm.size * n)
comm.cat("Original x:\n", quiet = TRUE)
comm.print(x, all.rank = TRUE)

x <- as.integer(x)
y <- spmd.alltoall.integer(x, integer(length(x)), n, n)
comm.cat("\nAlltoall y:\n", quiet = TRUE)
comm.print(y, all.rank = TRUE)

### Finish.
finalize()
