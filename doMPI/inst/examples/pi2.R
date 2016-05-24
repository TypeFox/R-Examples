# Here is yet another silly example that estimates pi.
# It demonstrates the use of the "timeout" function to limit
# the length of time to compute the estimate.  It is also
# yet another demonstration of using vector operations within
# a foreach loop.
#
# Note: This example will not work with the doMC or doSNOW
# parallel backends.  It only works with backends that submit
# tasks lazily, such as doMPI.  If used with "greedy" backends,
# like doMC or doSNOW, they could attempt to submit a huge number
# of tasks, resulting in some sort of nasty memory problems.

suppressMessages(library(doMPI))
suppressMessages(library(itertools))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Initialize variables
n <- 1000000  # length of vectors
t <- 60       # seconds to compute

# Create a "timeout" iterator that returns "n" for "t" seconds
timer <- timeout(irepeat(n), time=t)

# Define a ".final" function that calculates pi from estimates of pi/4
calc.pi <- function(x) {
  cat(sprintf('computed %d estimates of pi/4\n', length(x)))
  4 * mean(x)
}

# Compute pi in parallel
pi <- foreach(n=timer, .combine='c', .final=calc.pi) %dopar% {
  x <- runif(n=n, min=-0.5, max=0.5)
  y <- runif(n=n, min=-0.5, max=0.5)
  sum(sqrt(x*x + y*y) < 0.5) / n
}

# Display the result
cat(sprintf('Approximate value of pi after about %d seconds: %f\n', t, pi))

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
