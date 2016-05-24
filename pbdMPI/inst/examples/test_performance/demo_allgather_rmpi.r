library(Rmpi)
invisible(mpi.comm.dup(0, 1))

time.proc <- list()

time.proc$Robj <- system.time({
  for(i in 1:1000){
    y <- mpi.allgather.Robj(list(x = 1:10000))
  }
  mpi.barrier()
})

time.proc$matrix <- system.time({
  for(i in 1:1000){
    y <- mpi.allgather.Robj(matrix(1:10000, nrow = 100))
  }
  mpi.barrier()
})

if(mpi.comm.rank(1) == 0){
  print(time.proc)
}
mpi.quit()
