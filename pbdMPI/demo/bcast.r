### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
N <- 5
x <- (1:N) + N * .comm.rank
comm.cat("Original x:\n", quiet = TRUE)
comm.print(x, rank.print = 1)

y <- bcast(matrix(x, nrow = 1))
comm.cat("\nBcast matrix:\n", quiet = TRUE)
comm.print(y, rank.print = 1)

y <- bcast(as.integer(x))
comm.cat("\nBcast integer:\n", quiet = TRUE)
comm.print(y, rank.print = 1)

y <- bcast(as.double(x))
comm.cat("\nBcast double:\n", quiet = TRUE)
comm.print(y, rank.print = 1)

### Finish.
finalize()
