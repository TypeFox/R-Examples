### SHELL> mpiexec -np 2 Rscript --vanilla [...]_rmpi.r

library(Rmpi)
invisible(mpi.comm.dup(0, 1))

source("./01_setting")

x <- (1:N) + N * .comm.rank

time.proc <- list()

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- mpi.allreduce(as.integer(x), type = 1, op = "sum")
  }
  invisible(mpi.barrier())
})

time.proc$double <- system.time({ 
  for(i in 1:iter.total){
    y <- mpi.allreduce(as.double(x), type = 2, op = "sum")
  }
  invisible(mpi.barrier())
})

if(.comm.rank == 0){
  print(time.proc, quiet = TRUE)
}

mpi.quit()
