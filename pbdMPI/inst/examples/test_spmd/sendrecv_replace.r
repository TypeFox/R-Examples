### SHELL> mpiexec -np 2 Rscript --vanilla [...].r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

rank.dest <- as.integer((.comm.rank + 1) %% .comm.size)
rank.source <- as.integer((.comm.rank - 1) %% .comm.size)
comm.cat("-- Send to:\n", quiet = TRUE)
comm.print(rank.dest, all.rank = TRUE)

comm.cat("\n-- Receive from:\n", quiet = TRUE)
comm.print(rank.source, all.rank = TRUE)
comm.cat("\n", quiet = TRUE)

N <- 5
x <- (1:N) + N * .comm.rank

x <- sendrecv.replace(list(x))
comm.print(x, rank.print = 1)

x <- (1:N) + N * .comm.rank            ### since x has been replace by a list.
x <- sendrecv.replace(as.integer(x))
comm.print(x, rank.print = 1)

x <- sendrecv.replace(as.double(x))
comm.print(x, rank.print = 1)

finalize()
