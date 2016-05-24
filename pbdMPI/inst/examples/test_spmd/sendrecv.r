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
comm.cat("\n-- Result:\n", quiet = TRUE)

N <- 5
x <- (1:N) + N * .comm.rank

y <- sendrecv(list(x))
comm.print(y, rank.print = 1)

y <- sendrecv(as.integer(x), integer(N))
comm.print(y, rank.print = 1)

y <- sendrecv(as.double(x), double(N))
comm.print(y, rank.print = 1)

finalize()
