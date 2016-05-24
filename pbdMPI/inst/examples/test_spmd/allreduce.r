### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 5
x <- (1:N) + N * .comm.rank

y <- allreduce(matrix(as.integer(x), nrow = 1), op = "sum")
comm.print(y)
y <- allreduce(matrix(as.double(x), nrow = 1), op = "sum")
comm.print(y)

y <- allreduce(as.integer(x), op = "sum")
comm.print(y)
y <- allreduce(as.double(x), op = "sum")
comm.print(y)

y <- allreduce(as.integer(x), op = "prod")
comm.print(y)
y <- allreduce(as.double(x), op = "prod")
comm.print(y)

y <- allreduce(as.integer(x), op = "max")
comm.print(y)
y <- allreduce(as.double(x), op = "max")
comm.print(y)

y <- allreduce(as.integer(x), op = "min")
comm.print(y)
y <- allreduce(as.double(x), op = "min")
comm.print(y)

finalize()
