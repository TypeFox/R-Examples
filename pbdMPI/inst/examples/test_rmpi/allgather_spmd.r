### SHELL> mpiexec -np 2 Rscript --vanilla [...]_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

source("./01_setting")

x.total <- N * .comm.size
x <- (1:N) + N * .comm.rank

time.proc <- list()

time.proc$list <- system.time({
  for(i in 1:iter.total){
    y <- allgather(list(x))
  }
  barrier()
})

time.proc$matrix <- system.time({
  for(i in 1:iter.total){
    y <- allgather(matrix(x, nrow = sqrt(N)))
  }
  barrier()
})

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- allgather(as.integer(x), integer(x.total))
  }
  barrier()
})

time.proc$double <- system.time({ 
  for(i in 1:iter.total){
    y <- allgather(as.double(x), double(x.total))
  }
  barrier()
})

comm.print(time.proc, quiet = TRUE)

finalize()
