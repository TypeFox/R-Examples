### SHELL> mpiexec -np 2 Rscript --vanilla [...]_rmpi.r

library(Rmpi)
invisible(mpi.comm.dup(0, 1))

source("./01_setting")

x <- (1:N) + N * .comm.rank

time.proc <- list()

time.proc$Robj <- system.time({
  for(i in 1:iter.total){
    y <- mpi.bcast.Robj(matrix(x, nrow = sqrt(N)))
  }
  invisible(mpi.barrier())
})

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- mpi.bcast(as.integer(x), 1)
  }
  invisible(mpi.barrier())
})

time.proc$double <- system.time({
  for(i in 1:iter.total){
    y <- mpi.bcast(as.double(x), 2)
  }
  invisible(mpi.barrier())
})

if(.comm.rank == 0){
  print(time.proc)
}

mpi.quit()
