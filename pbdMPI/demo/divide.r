### SHELL> mpiexec -np 4 Rscript divide.r

### Initial.
suppressMessages(library(pbdMPI, quietly = TRUE))
init()
.comm.size <- comm.size()
.comm.rank <- comm.rank()

### Examples.
comm.cat(">>> block\n", quiet = TRUE)
jid <- get.jid(7, method = "block")
comm.print(jid, all.rank = TRUE)

comm.cat(">>> cycle\n", quiet = TRUE)
jid <- get.jid(7, method = "cycle")
comm.print(jid, all.rank = TRUE)

comm.cat(">>> block (all)\n", quiet = TRUE)
alljid <- get.jid(7, method = "block", all = TRUE)
comm.print(alljid)

comm.cat(">>> cycle (all)\n", quiet = TRUE)
alljid <- get.jid(7, method = "cycle", all = TRUE)
comm.print(alljid)

comm.cat(">>> block (all, reduced)\n", quiet = TRUE)
alljid.list <- get.jid(7, method = "block", all = TRUE, reduced = TRUE)
comm.print(alljid.list)

comm.cat(">>> cycle (all, reduced)\n", quiet = TRUE)
alljid.list <- get.jid(7, method = "cycle", all = TRUE, reduced = TRUE)
comm.print(alljid.list)


### Finish.
finalize()
