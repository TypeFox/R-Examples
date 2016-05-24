### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 1024
x <- (1:N) + N * .comm.rank

y <- reduce(as.integer(x), op = "sum")
comm.print(y[N - (1:5)])
y <- reduce(as.double(x), op = "sum")
comm.print(y[N - (1:5)])

y <- reduce(as.integer(x), op = "prod")
comm.print(y[N - (1:5)])
y <- reduce(as.double(x), op = "prod")
comm.print(y[N - (1:5)])

y <- reduce(as.integer(x), op = "max")
comm.print(y[N - (1:5)])
y <- reduce(as.double(x), op = "max")
comm.print(y[N - (1:5)])

y <- reduce(as.integer(x), op = "min")
comm.print(y[N - (1:5)])
y <- reduce(as.double(x), op = "min")
comm.print(y[N - (1:5)])

finalize()
