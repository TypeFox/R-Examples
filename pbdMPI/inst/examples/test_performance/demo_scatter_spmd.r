suppressMessages(library(pbdMPI, quietly = TRUE))
init()

time.proc <- list()

x <- split(1:10000, rep(1:2, 5000))
time.proc$Robj <- system.time({
  for(i in 1:1000){
    y <- scatter(x)
  }
  barrier()
})

time.proc$integer <- system.time({
  for(i in 1:1000){
    y <- scatter(as.integer(1:10000), integer(5000))
  }
  barrier()
})

comm.print(time.proc, quiet = TRUE)
finalize()
