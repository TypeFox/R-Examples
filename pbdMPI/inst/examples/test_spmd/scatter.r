### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 5
x.total <- N * .comm.size
x <- 1:x.total
y <- scatter(split(x, rep(1:.comm.size, N)))     ### return the element of list.
comm.print(y, rank.print = 1)
y <- scatter(as.integer(x), integer(N))
comm.print(y, rank.print = 1)
y <- scatter(as.double(x), double(N))
comm.print(y, rank.print = 1)

x.total <- (.comm.size + 1) * .comm.size / 2
x <- 1:x.total
x.count <- 1:.comm.size
y <- scatter(split(x, rep(x.count, x.count)))   ### return the element of list.
comm.print(y)
y <- scatter(as.integer(x), integer(.comm.rank + 1), as.integer(x.count))
comm.print(y, rank.print = 1)
y <- scatter(as.double(x), double(.comm.rank + 1), as.integer(x.count))
comm.print(y, rank.print = 1)

finalize()
