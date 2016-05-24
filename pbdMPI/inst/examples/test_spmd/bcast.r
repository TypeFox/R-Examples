### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 5
x <- (1:N) + N * .comm.rank
comm.print(x, rank.print = 1)

y <- bcast(matrix(x, nrow = 1))
comm.print(y, rank.print = 1)
y <- bcast(as.integer(x))
comm.print(y, rank.print = 1)
y <- bcast(as.double(x))
comm.print(y, rank.print = 1)

finalize()
