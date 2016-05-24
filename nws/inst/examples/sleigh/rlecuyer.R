library(nws)
library(rlecuyer)

# Worker function called via initRlecuyer on the master
initRlecuyerWorker <- function() {
  library(rlecuyer)

  # Fetch a "stream state" list
  slist <- nwsFetch(SleighNws, paste('RNGstream', SleighRank + 1, sep=''))

  # Recreate the random seed table as appropriate for this stream
  dim(slist$Cg) <- c(1, length(slist$Cg))
  dim(slist$Bg) <- c(1, length(slist$Bg))
  dim(slist$Ig) <- c(1, length(slist$Ig))
  AIP <- c(slist$Anti, slist$IncPrec)
  dim(AIP) <- c(1, 2)
  .lec.Random.seed.table <<- list(Cg=slist$Cg, Bg=slist$Bg, Ig=slist$Ig,
                                  AIP=AIP, name=slist$name)

  # Set the stream for the worker to use
  .lec.CurrentStream(slist$name)
  invisible(NULL)
}

# Initialize the sleigh to get random numbers using the rlecuyer package
initRlecuyer <- function(s, seed=rep(12345, 6)) {
  # Streams names, one for each sleigh worker
  snames <- paste('RNGstream', seq(length=workerCount(s)), sep='')

  # Call the init function so initRlecuyer can be called multiple times
  .lec.init()

  # Set the seed and create a stream for each worker
  .lec.SetPackageSeed(seed)
  .lec.CreateStream(snames)

  # Send the stream info needed by each worker, which will be
  # consumed by the worker function
  for (sn in snames)
    nwsStore(s@nws, sn, .lec.GetStateList(sn))
  eachWorker(s, initRlecuyerWorker)
  invisible(NULL)
}

# Simple test program for the initRlecuyer function
testRlecuyer <- function() {
  s <- sleigh(workerCount=10)
  seed <- 1:6
  initRlecuyer(s, seed)
  x <- unlist(eachWorker(s, runif, 10))

  # Make sure we get the same
  for (i in 1:10) {
    initRlecuyer(s, seed)
    y <- unlist(eachWorker(s, runif, 10))
    stopifnot(all(x == y))
  }

  stopSleigh(s)
  x
}
