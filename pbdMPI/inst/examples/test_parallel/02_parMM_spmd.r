### This example is analog to "parMM_par.r", and one can run it by the command
### SHELL> mpiexec -np 2 Rscript --vanilla 02_parMM_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

set.seed(123)
x <- matrix(rnorm(1000000), 1000)

parMM.spmd <- function(x, y){
  id <- get.jid(nrow(x))
  do.call(rbind, allgather(x[id,] %*% y))
}
time.proc <- system.time(replicate(10, parMM.spmd(x, x)))
comm.print(time.proc)

finalize()
