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
comm.print(x, all.rank = TRUE)

y <- allreduce(matrix(x, nrow = 1), op = "sum")
comm.cat("\nAllreduce sum:\n", quiet = TRUE)
comm.print(y)

y <- allreduce(x, op = "prod")
comm.cat("\nAllreduce prod:\n", quiet = TRUE)
comm.print(y)

y <- allreduce(x, op = "max")
comm.cat("\nAllreduce max:\n", quiet = TRUE)
comm.print(y)

y <- allreduce(x, op = "min")
comm.cat("\nAllreduce min:\n", quiet = TRUE)
comm.print(y)

### Finish.
finalize()
