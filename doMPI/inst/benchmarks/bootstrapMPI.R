# NOTE:  I consider this to be a useful benchmark for comparing
# different parallel systems, but it is important to keep in
# mind that the tasks are rather small.  That is useful for
# finding out how much overhead different parallel programming
# systems add to task execution.  But because the tasks are small,
# it isn't really a great candidate for parallel programming.
# However, it isn't ridiculous either, and would be a very good
# example using hardware from a decade ago.
#
# Also note that because this is a benchmark, I am not actually
# combining the task results, but am throwing them away.  That
# allows me to test just the speed of the parallel programming
# system, not the speed of the combine mechanism in the foreach
# package.
#
# TODO Document the usage

suppressMessages(library(doMPI))
suppressMessages(library(getopt))

main <- function(args) {
  # Define the command line options
  spec <- matrix(c('verbose',    'v', '0', 'logical',
                   'profile',    'p', '0', 'logical',
                   'emulate',    'e', '0', 'logical',
                   'force',      'f', '0', 'logical',
                   'threshold',  't', '1', 'integer',
                   'cores',      'c', '1', 'integer',
                   'nochunking', 'k', '0', 'logical',
                   'sequential', 's', '0', 'logical'),
                 ncol=4, byrow=TRUE)

  # Process the command line
  options <- getopt(spec, opt=args, command='bootstrapMPI.R')

  # Define the default value of the options
  opt <- list(verbose=FALSE, profile=FALSE, emulate=FALSE, force=FALSE,
              threshold=800, cores=1, nochunking=FALSE, sequential=FALSE)

  # Merge the command line with the default values
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
    cat(sprintf("Chunking is %s\n", if (opt$nochunking) "off" else "on"))

    x <- iris[which(iris[,5] != "setosa"), c(1,5)]
    trials <- 10000

    chunkSize <- if (opt$nochunking) 1 else ceiling(trials / getDoParWorkers())
    opts <- list(chunkSize=chunkSize, profile=opt$profile,
                 forcePiggyback=opt$force,
                 bcastThreshold=opt$threshold)
    trash <- function(...) NULL

    ptime <- system.time({
      foreach(icount(trials), .combine=trash, .multicombine=TRUE,
                   .options.mpi=opts) %dopar% {
        ind <- sample(100, 100, replace=TRUE)
        result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
        coefficients(result1)
      }
    })[3]

    cat(sprintf('Parallel time using doMPI on %d workers: %f\n',
                getDoParWorkers(), ptime))

    # We must shutdown the cluster before running the sequential benchmark
    # in case any cluster workers are running on the local machine
    closeCluster(cl)

    if (opt$sequential) {
      stime <- system.time({
        foreach(icount(trials), .combine=trash, .multicombine=TRUE) %do% {
          ind <- sample(100, 100, replace=TRUE)
          result1 <- glm(x[ind,2]~x[ind,1], family=binomial(logit))
          coefficients(result1)
        }
      })[3]

      cat(sprintf('Sequential time: %f\n', stime))
      cat(sprintf('Speed up for %d workers: %f\n',
                  getDoParWorkers(), round(stime / ptime, digits=2)))
    }
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
