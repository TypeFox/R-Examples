suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Initialize variables
n <- 10000
opts <- list(chunkSize=100)

# Perform simulations in parallel
r <- foreach(1:n, .combine='c', .options.mpi=opts) %dopar% {
  x <- ts(arima.sim(list(order=c(1,0,0), ar=-0.9), n=360), start=1975, freq=12)
  y <- aggregate(x, nfreq=4, sum)
  arima(y, order=c(1,0,0))$model$phi
}

# Print a summary of the resulting vector
print(summary(r))

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
