suppressPackageStartupMessages(library(pbdMPI, quietly=TRUE))
init()

sleeptime <- 1

time <- system.time({
  if (comm.rank()==0)
    Sys.sleep(sleeptime)
  
  barrier()
})[3]

test <- abs(time - sleeptime) <= .1 # sleep isn't really guaranteed
test <- comm.all(test)
comm.print(test)

finalize()
