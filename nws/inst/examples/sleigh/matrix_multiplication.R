# matrix multiplication
# each worker computes one block column of the result matrix

library(nws)

as.charlist <- function(s) {
  x <- strsplit(s, " +", extended=TRUE)[[1]]
  x[x != '']
}

as.intlist <- function(s) {
  x <- as.integer(strsplit(s, " +", extended=TRUE)[[1]])
  x[!is.na(x)]
}

######################## Configuration Section #########################
NODELIST <- as.charlist(Sys.getenv('NODELIST', ''))
NUMWORKERS <- if (length(NODELIST) > 0) length(NODELIST) else 4
NUMWORKERS <- as.intlist(Sys.getenv('NUMWORKERS', NUMWORKERS))
SEQUENTIAL <- as.logical(Sys.getenv('SEQUENTIAL', 'FALSE'))
VERBOSE <- as.logical(Sys.getenv('VERBOSE', 'FALSE'))
NUMROWS <- as.integer(Sys.getenv('NUMROWS', 1200))
########################################################################

if (is.na(SEQUENTIAL)) {
  stop('bad value for SEQUENTIAL: ', Sys.getenv('SEQUENTIAL'))
}

# worker initialization function
initWorker <- function() {
  assign('a', nwsFind(SleighUserNws, 'a'), envir=globalenv())
  NULL
}

# worker task function
mvmult <- function(y) {
  a %*% y
}

# worker clean up function
cleanWorker <- function() {
  rm('a', envir=globalenv())
}

# block column extract function
colextract <- function(j, n, x) {
  '['(x, , seq(j, length.out=min(n, dim(x)[2] - j + 1)), drop=FALSE)
}

# parallel matrix multiplication function
matMult <- function(s, a, b) {
  cat('Storing matrix\n')
  nwsStore(s@userNws, 'a', a)
  cat('Finding matrix\n')
  eachWorker(s, initWorker)
  cat('Deleting matrix variable\n')
  nwsDeleteVar(s@userNws, 'a')

  cat('Calling mvmult on workers\n')
  blocksize <- ceiling(dim(b)[2] / workerCount(s))
  lb <- lapply(seq(1, dim(b)[2], by=blocksize), colextract, blocksize, b)
  p.time <- system.time(result <- eachElem(s, mvmult, list(lb)))
  cat('eachElem time: ')
  print(p.time)
  cat('Deleting global variable on workers\n')
  eachWorker(s, cleanWorker)
  cat('Constructing result matrix\n')
  do.call('cbind', result)
}

# generate two random matrices
n <- NUMROWS * NUMROWS
A <- matrix(rnorm(n), NUMROWS, NUMROWS)
B <- matrix(rnorm(n), NUMROWS, NUMROWS)

if (SEQUENTIAL) {
  # multiply the matrices sequentially
  cat('Running sequentially\n')
  seq.time <- system.time(D <- A %*% B)
  cat('Sequential time:\n')
  print(seq.time)
}

# multiply the matrices in parallel
for (numWorkers in NUMWORKERS) {
  cat('Creating a sleigh\n')
  s <- if (length(NODELIST) > 0)
         sleigh(nodeList=rep(NODELIST, length.out=numWorkers),
                launch=sshcmd, verbose=VERBOSE)
       else
         sleigh(workerCount=numWorkers, verbose=VERBOSE)
  cat('Waiting for the sleigh to start\n')
  stat <- status(s, TRUE, 60)
  cat('Running on', stat$numWorkers, 'workers\n')
  par.time <- system.time(C <- matMult(s, A, B))
  cat('Parallel time for', stat$numWorkers, 'workers:\n')
  print(par.time)
  if (SEQUENTIAL) {
    cat('Result comparison: ')
    print(all.equal(C, D))
  }
  rm(C)
  cat('Shutting down the sleigh\n')
  stopSleigh(s)
}
