suppressMessages(library(pbdMPI, quietly = TRUE))
init()

time.proc <- list()

time.proc$default <- system.time({
  for(i in 1:1000){
    y <- allgather(list(x = 1:10000))
  }
  barrier()
})

time.proc$matrix <- system.time({
  for(i in 1:1000){
    y <- allgather(matrix(1:10000, nrow = 100))
  }
  barrier()
})

comm.print(time.proc, quiet = TRUE)
finalize()
