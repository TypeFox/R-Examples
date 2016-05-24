### SHELL> mpiexec -np 2 Rscript --vanilla [...]_rmpi.r

library(Rmpi)
invisible(mpi.comm.dup(0, 1))

source("./01_setting")

x <- 1:(N * .comm.size)
x.list <- split(x, rep(1:.comm.size, N))
x.count <- rep(N, .comm.size)

time.proc <- list()

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- mpi.scatterv(as.integer(x), as.integer(x.count), 1, integer(N))
  }
  invisible(mpi.barrier())
})

time.proc$double <- system.time({
  for(i in 1:iter.total){
    y <- mpi.scatterv(as.double(x), as.integer(x.count), 2, double(N))
  }
  invisible(mpi.barrier())
})

if(.comm.rank == 0){
  print(time.proc)
}

mpi.quit()
