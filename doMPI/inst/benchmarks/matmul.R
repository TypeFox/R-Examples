# This isn't an example.  It's a benchmark, so it has a number
# of options to control the execution of the benchmark.
# That makes it too complicated for an example, although the
# complexity has nothing to do with parallel computing.
#
# TODO Document the usage

suppressMessages(library(doMPI))
suppressMessages(library(getopt))

# Define a parallel matrix multiply function
matmul <- function(x, y, profile=FALSE, forcePiggyback=FALSE,
                   bcastThreshold=800) {
  n <- ceiling(ncol(y) / getDoParWorkers())
  opts <- list(profile=profile, forcePiggyback=forcePiggyback,
               bcastThreshold=bcastThreshold)
  foreach(yc=iter(y, by='column', chunksize=n), .combine='cbind',
          .options.mpi=opts) %dopar% {
    x %*% yc
  }
}

# Main program executed by master and workers
main <- function(args) {
  spec <- matrix(c('verbose',   'v', '0', 'logical',
                   'profile',   'p', '0', 'logical',
                   'emulate',   'e', '0', 'logical',
                   'force',     'f', '0', 'logical',
                   'threshold', 't', '1', 'integer',
                   'cores',     'c', '1', 'integer'), ncol=4, byrow=TRUE)
  options <- getopt(spec, opt=args, command='matmul.R')
  opt <- list(verbose=FALSE, profile=FALSE, emulate=FALSE,
              force=FALSE, threshold=800, cores=1)
  opt[names(options)] <- options

  # Check if this is the master or a cluster worker
  if (mpi.comm.rank(0) > 0) {
    # This is a cluster worker
    wfile <- sprintf("MPI_%d_%s.log", mpi.comm.rank(0), Sys.info()[['user']])
    outfile <- if (opt$verbose) wfile else "/dev/null"
    sinkWorkerOutput(outfile)
    cl <- openMPIcluster(bcast=!opt$emulate)
    dompiWorkerLoop(cl, cores=opt$cores, verbose=opt$verbose)
  } else {
    # Create and register an MPI cluster
    cl <- startMPIcluster(bcast=!opt$emulate, verbose=opt$verbose)
    registerDoMPI(cl)

    # Display a summary of how we're going to run the benchmark
    if (opt$force) {
      cat("Broadcasting is disabled: job data will always be piggy-backed\n")
    } else {
      if (opt$emulate) {
        cat("Using emulated broadcast for job data\n")
      } else {
        cat("Using true MPI broadcast for job data\n")
      }
      cat(sprintf("Piggy-back/broadcast threshold is %d\n", opt$threshold))
    }

    # Create some matrices
    m <- 2000; n <- 1200; p <- 1200
    x <- matrix(rnorm(m * n), m, n)
    y <- matrix(rnorm(n * p), n, p)

    # Execute matmul and report the time
    stime <- proc.time()[3]
    z <- matmul(x, y, opt$profile, opt$force, opt$threshold)
    etime <- proc.time()[3] - stime
    cat(sprintf('Time for matrix multiply with %s and %d workers: %f\n',
                getDoParName(), getDoParWorkers(), etime))

    # Shutdown the cluster
    closeCluster(cl)
  }
}

# Catch errors so we can call mpi.quit
tryCatch({
  main(commandArgs(trailingOnly=TRUE))
},
error=function(e) {
  cat(sprintf('%s\n', conditionMessage(e)))
})

mpi.quit()
