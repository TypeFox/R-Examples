# This example shows how to use the .combine option to write results
# to a file created by the master process.  It demonstrates how the
# foreach .init and .final arguments are used as well, and how the
# first argument to the .combine function can be different than the
# subsequent arguments.
#
# See initEnvir.R for a similar example, but which has the
# workers write the data themselves.

suppressMessages(library(doMPI))

# Create and register an MPI cluster
cl <- startMPIcluster()
registerDoMPI(cl)

# Define a .combine function that writes its arguments to a file
writeResults <- function(fobj, ...) {
  args <- list(...)
  cat(sprintf('writeResults called with %d results\n', length(args)))
  d <- do.call('rbind', args)
  d$y <- signif(d$y, 5)
  write.table(d, file=fobj, col.names=FALSE, row.names=FALSE, sep=",",
              quote=FALSE, eol="\r\n", na=" ", dec=".")
  flush(fobj)  # useful for debugging, could be slow
  fobj  # the return value must be the file object
}

# Driver function
main <- function(fname) {
  # Create the output file that we will specify via the .init argument
  f <- file(fname, 'w')
  dummy <- data.frame(i=integer(0), x=integer(0), y=double(0))
  write.table(dummy, file=f, col.names=TRUE, row.names=FALSE, sep=",",
              quote=FALSE, eol="\r\n", na=" ", dec=".")
  flush(f)

  foreach(i=icount(100), .init=f, .final=close, .combine=writeResults,
          .maxcombine=20) %dopar% {
    data.frame(i=i, x=1:4, y=rnorm(4))
  }
}

# Call the driver function
main('dompi-1.csv')

# Shutdown the cluster and quit
closeCluster(cl)
mpi.quit()
