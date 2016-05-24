### This example is analog to "glm_par.r", and one can run it by the command
### SHELL> mpiexec -np 2 Rscript --vanilla 01_mclapply_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

time.proc <- system.time({
  id <- get.jid(32)
  ret <- unlist(lapply(id, function(i) sum(rnorm(1e7))))
  allgather(ret, unlist = TRUE)
})
comm.print(time.proc)

finalize()
