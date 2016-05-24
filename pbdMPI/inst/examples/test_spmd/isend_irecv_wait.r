### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 5
x <- (1:N) + N * .comm.rank

if(.comm.rank == 0){
  isend(list(x))
}
if(.comm.rank == 1){
  y <- irecv(list(x))
}
wait()
comm.print(y, rank.print = 1)

if(.comm.rank == 0){
  isend(as.integer(x))
}
if(.comm.rank == 1){
  y <- irecv(as.integer(x))
}
z <- waitany(as.integer(comm.size()))
comm.print(z, rank.print = 1)
comm.print(y, rank.print = 1)

if(.comm.rank == 0){
  isend(as.double(x))
}
if(.comm.rank == 1){
  y <- irecv(as.double(x))
}
z <- waitsome(as.integer(comm.size()))
comm.print(z, rank.print = 1)
comm.print(y, rank.print = 1)

waitall(as.integer(comm.size()))

finalize()
