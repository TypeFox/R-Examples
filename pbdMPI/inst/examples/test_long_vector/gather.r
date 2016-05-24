### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 1024
x.total <- N * .comm.size
x <- (1:N) + N * .comm.rank
y <- gather(as.integer(x), integer(x.total))
comm.print(y[N - (1:5)])
y <- gather(as.double(x), double(x.total))
comm.print(y[N - (1:5)])

finalize()
