# This example is an R version of an example in "Getting Started with MATLAB",
# on "Graphing the sinc function".  It's an interesting example of using
# recycling in vector operations to write a clean, parallel version.
# And it also makes a pretty picture, which is why there is also a demo
# version of it included with doMPI.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Define the grid over which to compute the sinc function
x <- seq(-20, 20, by=0.25)

# Compute the sinc function in parallel
v <- foreach(y=x, .combine="cbind") %dopar% {
  r <- sqrt(x^2 + y^2) + .Machine$double.eps
  sin(r) / r
}

# Display the results with a perspective plot
persp(x, x, v)

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
