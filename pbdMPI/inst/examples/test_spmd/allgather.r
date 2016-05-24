### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 5
x.total <- N * .comm.size
x <- (1:N) + N * .comm.rank
y <- allgather(matrix(x, nrow = 1))
comm.print(y)
y <- allgather(as.integer(x), integer(x.total))
comm.print(y)
y <- allgather(as.double(x), double(x.total))
comm.print(y)

x <- 1:(.comm.rank + 1)
x.total <- (.comm.size + 1) * .comm.size / 2
x.count <- 1:.comm.size
y <- allgather(matrix(x, nrow = 1))
comm.print(y)
y <- allgather(as.integer(x), integer(x.total), as.integer(x.count))
comm.print(y)
y <- allgather(as.double(x), double(x.total), as.integer(x.count))
comm.print(y)

finalize()
