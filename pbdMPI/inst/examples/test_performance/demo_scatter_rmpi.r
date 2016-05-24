library(Rmpi)
invisible(mpi.comm.dup(0, 1))

time.proc <- list()

x <- split(1:10000, rep(1:2, 5000))
time.proc$Robj <- system.time({
  for(i in 1:1000){
    y <- mpi.scatter.Robj(x)
  }
  invisible(mpi.barrier())
})

time.proc$integer <- system.time({
  for(i in 1:1000){
    y <- mpi.scatter(as.integer(1:10000), 1, integer(5000))
  }
  invisible(mpi.barrier())
})

if(mpi.comm.rank(1) == 0){
  print(time.proc)
}
mpi.quit()
