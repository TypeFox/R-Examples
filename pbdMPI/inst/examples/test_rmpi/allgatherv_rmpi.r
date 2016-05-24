### SHELL> mpiexec -np 2 Rscript --vanilla [...]_rmpi.r

library(Rmpi)
invisible(mpi.comm.dup(0, 1))

source("./01_setting")

x.total <- N * .comm.size
x <- (1:N) + N * .comm.rank
x.count <- rep(N, .comm.size)

time.proc <- list()

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- mpi.allgatherv(as.integer(x), 1, integer(x.total), as.integer(x.count))
  }
  invisible(mpi.barrier())
})

time.proc$double <- system.time({
  for(i in 1:iter.total){
    y <- mpi.allgatherv(as.double(x), 2, double(x.total), as.integer(x.count))
  }
  invisible(mpi.barrier())
})

if(.comm.rank == 0){
  print(time.proc)
}

mpi.quit()
