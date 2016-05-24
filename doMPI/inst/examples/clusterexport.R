# This example shows how to export data to the workers so that it
# can be reused across multiple foreach loops using the standard
# 'attach' function and the doMPI-specific 'initEnvir' option.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Initialize variables
x <- 7

# Export 'x' to all of the workers for use across multiple foreach loops
initEnvir <- function(envir, e) attach(e, name='dompi_example')
e <- new.env(parent=emptyenv())
e$x <- x
mpiopts <- list(initEnvir=initEnvir, initArgs=list(e=e))
ignore <- foreach(icount(getDoParWorkers()), .options.mpi=mpiopts) %dopar% NULL

# Use the 'x' that is attached
r <- foreach(i=1:10, .combine='c', .noexport='x') %dopar% {
  i * x
}
print(r)

# Modify 'x'
x <- 6

# Use the new value of 'x' which is auto-exported
r <- foreach(i=1:10, .combine='c') %dopar% {
  i * x
}
print(r)

# Use the original 'x' that is still attached in all workers
r <- foreach(i=1:10, .combine='c', .noexport='x') %dopar% {
  i * x
}
print(r)

# Unexport 'x'
initEnvir <- function(envir, e) detach(name='dompi_example')
mpiopts <- list(initEnvir=initEnvir)
ignore <- foreach(icount(getDoParWorkers()), .options.mpi=mpiopts) %dopar% NULL

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
