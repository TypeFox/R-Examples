### SHELL> mpiexec -np 2 Rscript --vanilla [...]_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

source("./01_setting")

x <- (1:N) + N * .comm.rank

time.proc <- list()

time.proc$array <- system.time({
  for(i in 1:iter.total){
    y <- allreduce(matrix(as.integer(x), nrow = sqrt(N)), op = "sum")
  }
  barrier()
})

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- allreduce(as.integer(x), integer(N), op = "sum")
  }
  barrier()
})

time.proc$double <- system.time({ 
  for(i in 1:iter.total){
    y <- allreduce(as.double(x), double(N), op = "sum")
  }
  barrier()
})

time.proc$integer.nobuffer <- system.time({
  for(i in 1:iter.total){
    y <- allreduce(as.integer(x), op = "sum")
  }
  barrier()
})

time.proc$double.nobuffer <- system.time({ 
  for(i in 1:iter.total){
    y <- allreduce(as.double(x), op = "sum")
  }
  barrier()
})

comm.print(time.proc, quiet = TRUE)

finalize()
