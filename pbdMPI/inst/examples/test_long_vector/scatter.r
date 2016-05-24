### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 1024
x.total <- N * .comm.size
x <- 1:x.total
y <- scatter(as.integer(x), integer(N))
comm.print(y[N - (1:5)], rank.print = 1)
y <- scatter(as.double(x), double(N))
comm.print(y[N - (1:5)], rank.print = 1)

finalize()
