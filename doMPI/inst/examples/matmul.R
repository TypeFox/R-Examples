# This example shows a very simple way of performing a simple,
# parallel matrix multiply.  Of course, it's very inefficient,
# since matrix multiplication is so fast on modern computers.
# But it's something of a traditional in the parallel computing
# world, and it shows how to use the "iter" method on a matrix
# to iterate over block columns of the matrix "y".

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Define a parallel matrix multiply function
matmul <- function(x, y) {
  # one task per worker for maximum granularity, but it's still hopeless
  n <- ceiling(ncol(y) / getDoParWorkers())
  foreach(yc=iter(y, by='column', chunksize=n), .combine='cbind') %dopar% {
    x %*% yc
  }
}

# Create some matrices
m <- 6; n <- 5; p <- 4
x <- matrix(rnorm(m * n), m, n)
y <- matrix(rnorm(n * p), n, p)

# Execute matmul and then display and check the result
z <- matmul(x, y)
print(z)
print(identical(z, x %*% y))

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
