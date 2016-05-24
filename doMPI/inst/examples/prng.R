# This example shows how to setup the cluster workers to use parallel
# random number generation using the new setRngDoMPI function.  It
# causes each worker to produce independent random numbers, but the
# results aren't necessarily reproducible since the tasks can be
# executed by different workers on different runs.  The results will
# also be different for different numbers of workers.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

trials <- 4
w <- clusterSize(cl)

cat(sprintf("Only first %d results are guaranteed repeatable\n", w))
fun <- function(trial, comm) {
  # Initialize parallel RNG
  setRngDoMPI(cl, seed=42)

  foreach(sleep=irunif(1, max=5, count=5*w),
          .combine='rbind') %dopar% {
    Sys.sleep(sleep)  # Randomize task length
    data.frame(rank=mpi.comm.rank(comm), result=as.integer(runif(1, max=1000)))
  }
}
r <- lapply(1:trials, fun, cl$comm)
print(do.call('cbind', r))

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
