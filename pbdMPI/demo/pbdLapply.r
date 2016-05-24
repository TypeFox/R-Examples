### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
N <- 100
x <- split((1:N) + N * .comm.rank, rep(1:10, each = 10))
y <- pbdLapply(x, sum, pbd.mode = "spmd")
comm.print(unlist(y), all.rank = TRUE)

y <- pbdLapply(x, sum)
comm.print(unlist(y), all.rank = TRUE)

y <- pbdLapply(x, sum, bcast = TRUE)
comm.print(unlist(y), all.rank = TRUE)

### Finish.
finalize()
