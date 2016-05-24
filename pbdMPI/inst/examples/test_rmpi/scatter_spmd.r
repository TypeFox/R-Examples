### SHELL> mpiexec -np 2 Rscript --vanilla [...]_spmd.r

suppressMessages(library(pbdMPI, quietly = TRUE))
init()

source("./01_setting")

x <- 1:(N * .comm.size)
x.list <- split(x, rep(1:.comm.size, N))

time.proc <- list()

time.proc$default <- system.time({
  for(i in 1:iter.total){
    y <- scatter(x.list)
  }
  barrier()
})

time.proc$integer <- system.time({
  for(i in 1:iter.total){
    y <- scatter(as.integer(x), integer(N))
  }
  barrier()
})

time.proc$double <- system.time({
  for(i in 1:iter.total){
    y <- scatter(as.double(x), double(N))
  }
  barrier()
})

comm.print(time.proc, quiet = TRUE)

finalize()
