# This example shows how to export data to the workers so that it
# can be reused across multiple foreach loops by taking advantage
# of the mpirun.  This is particularly useful when you want to
# load a large amount of data once to initialize the workers.

suppressMessages(library(doMPI))

# Initialize variables to be "exported" before calling 'startMPIcluster'.
# Note that you can also use the 'attach' function rather than
# global variables.
if (TRUE) {
  x <- 7
} else {
  e <- new.env()
  e$x <- 7
  attach(e, name='exported_data')
  rm(e)
}

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Verify that this script was executed correctly
if (mpi.comm.size(0) == 1) {
  cat('\nYou must start the workers via mpirun with a command such as:\n')
  cat('  mpirun -n 4 R --slave -f clusterexport2.R\n\n')
  cat('This technique for initializing workers depends on being\n')
  cat('executed in "Single Program/Multiple Data" style.\n')
} else {
  # Use the 'x' in .GlobalEnv
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

  # Use the original 'x' that is still in .GlobalEnv
  r <- foreach(i=1:10, .combine='c', .noexport='x') %dopar% {
    i * x
  }
  print(r)
}

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
