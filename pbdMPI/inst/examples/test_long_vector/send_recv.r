### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

N <- 1024
x <- (1:N) + N * .comm.rank

if(.comm.rank == 0){
  send(as.integer(x))
} else if(.comm.rank == 1){
  y <- recv(as.integer(x))
}
comm.print(y[N - (1:5)], rank.print = 1)

if(.comm.rank == 0){
  send(as.double(x))
} else if(.comm.rank == 1){
  y <- recv(as.double(x))
}
comm.print(y[N - (1:5)], rank.print = 1)

finalize()
