### SHELL> mpiexec -np 2 Rscript --vanilla [...]_rmpi.r

library(Rmpi)
invisible(mpi.comm.dup(0, 1))

source("./01_setting")

x <- (1:N) + N * .comm.rank

time.proc <- list()

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- mpi.reduce(as.integer(x), 1, op = "sum")
  }
  invisible(mpi.barrier())
})

time.proc$double <- system.time({ 
  for(i in 1:iter.total){
    y <- mpi.reduce(as.double(x), 2, op = "sum")
  }
  invisible(mpi.barrier())
})

if(.comm.rank == 0){
  print(time.proc, quiet = TRUE)
}

mpi.quit()
