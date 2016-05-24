### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
N <- 5
x <- 1:(.comm.rank + 1)
x.total <- (.comm.size + 1) * .comm.size / 2
x.count <- 1:.comm.size
comm.cat("Original x:\n", quiet = TRUE)
comm.print(x, all.rank = TRUE)

y <- gather(matrix(x, nrow = 1))
comm.cat("\nGather matrix:\n", quiet = TRUE)
comm.print(y)

y <- gather(as.integer(x), integer(x.total), as.integer(x.count))
comm.cat("\nGatherv integer:\n", quiet = TRUE)
comm.print(y)

y <- gather(as.double(x), double(x.total), as.integer(x.count))
comm.cat("\nGatherv double:\n", quiet = TRUE)
comm.print(y)

comm.cat("\n-- Print from rank 1:\n", rank.print = 1, quiet = TRUE)
comm.print(y, rank.print = 1)

### Finish.
finalize()
