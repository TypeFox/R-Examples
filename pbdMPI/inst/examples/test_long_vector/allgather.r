### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 1024
x.total <- N * .comm.size
x <- (1:N) + N * .comm.rank
y <- allgather(as.integer(x), integer(x.total))
comm.print(y[N - (1:5)])
y <- allgather(as.double(x), double(x.total))
comm.print(y[N - (1:5)])

x <- 1:(.comm.rank + 1)
x.total <- (.comm.size + 1) * .comm.size / 2
x.count <- 1:.comm.size
y <- allgather(as.integer(x), integer(x.total), as.integer(x.count))
comm.print(y[N - (1:5)])
y <- allgather(as.double(x), double(x.total), as.integer(x.count))
comm.print(y[N - (1:5)])

finalize()
