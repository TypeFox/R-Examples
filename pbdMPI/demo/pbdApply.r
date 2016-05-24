### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
N <- 100
x <- matrix((1:N) + N * .comm.rank, ncol = 10)
comm.print(x)

y <- pbdApply(x, 1, sum, pbd.mode = "mw")
comm.print(y)

y <- pbdApply(x, 1, sum, pbd.mode = "spmd")
comm.print(y)

### Finish.
finalize()
