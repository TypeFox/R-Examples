### SHELL> mpiexec -np 2 Rscript --vanilla [...]_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

source("./01_setting")

x.total <- N * .comm.size
x <- (1:N) + N * .comm.rank
x.count <- rep(N, .comm.size)

time.proc <- list()

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- allgather(as.integer(x), integer(x.total), as.integer(x.count))
  }
  barrier()
})

time.proc$double <- system.time({ 
  for(i in 1:iter.total){
    y <- allgather(as.double(x), double(x.total), as.integer(x.count))
  }
  barrier()
})

comm.print(time.proc, quiet = TRUE)

finalize()
