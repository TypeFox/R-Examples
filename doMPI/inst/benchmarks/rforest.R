# This isn't an example.  It's a benchmark, so it has a number
# of options to control the execution of the benchmark.
# That makes it too complicated for an example, although the
# complexity has nothing to do with parallel computing.
#
# TODO Document the usage

suppressMessages(library(doMPI))
suppressMessages(library(randomForest))
suppressMessages(library(getopt))

# Define a parallel randomForest function
rforest <- function(x, y=NULL, xtest=NULL, ytest=NULL, ntree=500, ...,
                    profile=FALSE, forcePiggyback=FALSE, bcastThreshold=800) {
  comm <- getDoMpiCluster()$comm
  if (is.null(comm)) {
    stop('error getting communicator')
  }
  winit <- function(envir, dimen, comm) {
    library(randomForest)
    envir$x <- mpi.bcast(double(dimen[1] * dimen[2]), 2, rank=0, comm=comm)
    dim(envir$x) <- dimen
  }
  wargs=list(dimen=dim(x), comm=comm)
  minit <- function() {
    mpi.bcast(x, 2, rank=0, comm=comm)
  }
  opts <- list(profile=profile, forcePiggyback=forcePiggyback,
               bcastThreshold=bcastThreshold,
               initEnvir=winit, initArgs=wargs,
               initEnvirMaster=minit)

  foreach(i=idiv(ntree, chunks=getDoParWorkers()),
          .combine='combine', .multicombine=TRUE, .inorder=FALSE,
          .noexport='x', .options.mpi=opts) %dopar% {
    randomForest:::randomForest.default(x, y, xtest, ytest, ntree=i, ...)
  }
}

# Main program executed by master and workers
main <- function(args) {
  spec <- matrix(c('verbose',   'v', '0', 'logical',
                   'profile',   'p', '0', 'logical',
                   'emulate',   'e', '0', 'logical',
                   'force',     'f', '0', 'logical',
                   'threshold', 't', '1', 'integer',
                   'cores',     'c', '1', 'integer',
                   'rows',      'm', '1', 'integer',
                   'cols',      'n', '1', 'integer',
                   'ntree',     'r', '1', 'integer',
                   'importance','i', '0', 'logical'), ncol=4, byrow=TRUE)
  options <- getopt(spec, opt=args, command='rforest.R')
  opt <- list(verbose=FALSE, profile=FALSE, emulate=FALSE,
              force=FALSE, threshold=800, cores=1, rows=200, cols=100,
              ntree=500, importance=FALSE)
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

    # Create a matrix and factor as input
    m <- opt$rows; n <- opt$cols
    x <- matrix(rnorm(m * n), m, n)
    y <- gl(10, m/10)

    # Execute rforest and report the time
    stime <- proc.time()[3]
    rfit <- rforest(x, y, ntree=opt$ntree, importance=opt$importance,
                    profile=opt$profile, forcePiggyback=opt$force,
                    bcastThreshold=opt$threshold)
    print(rfit)
    etime <- proc.time()[3] - stime
    cat(sprintf('Time for randomForest with %s and %d workers: %f\n',
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
