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
NUMTASKS <- as.integer(Sys.getenv('NUMTASKS', 1000))
NUMSAMPLES <- as.integer(Sys.getenv('NUMSAMPLES', 20000))
NUMSTOCKS <- as.integer(Sys.getenv('NUMSTOCKS', 25))
CHUNKSIZE <- as.integer(Sys.getenv('CHUNKSIZE', 10))
LOADFACTOR <- as.integer(Sys.getenv('LOADFACTOR', 3))
########################################################################

if (is.na(SEQUENTIAL)) {
  stop('bad value for SEQUENTIAL: ', Sys.getenv('SEQUENTIAL'))
}

# hack to make the computation function work sequentially
SleighRank <- 1

# randomly generate the mean and sd that describe each stock
smean <- rnorm(NUMSTOCKS, mean=10.0, sd=1.0)
ssd <- rnorm(NUMSTOCKS, mean=3.0, sd=0.5)
stocks <- data.frame(mean=smean, sd=ssd)

# this is the task function, called via eachElem and lapply
fun <- function(numSamples, numStocks) {
  # generate the weights vector
  t <- runif(numStocks)
  w <- t / sum(t)

  # generate random stock returns matrix
  rnormWrapper <- function(i)
    rnorm(numSamples, mean=stocks$mean[[i]], sd=stocks$sd[[i]])
  s <- do.call(rbind, lapply(1:numStocks, rnormWrapper))

  # do the computation and return the results
  r <- drop(w %*% s)
  c(mean(r), var(r), SleighRank)
}

# do the work sequentially
if (SEQUENTIAL) {
  cat('Running sequentially\n')
  seq.time <- system.time(
    r2 <- lapply(rep(NUMSAMPLES, NUMTASKS), fun, NUMSTOCKS))
  cat('Sequential time:\n')
  print(seq.time)
}

# do the work in parallel
for (numWorkers in NUMWORKERS) {
  cat('Creating a sleigh\n')
  s <- if (length(NODELIST) > 0)
         sleigh(nodeList=rep(NODELIST, length.out=numWorkers),
                launch=sshcmd, verbose=VERBOSE)
       else
         sleigh(workerCount=numWorkers, verbose=VERBOSE)
  cat('Waiting for the sleigh to start\n')
  stat <- status(s, TRUE, 20)
  cat('Running on', stat$numWorkers, 'workers\n')
  tmp <- eachWorker(s, function(g1) {stocks <<- g1; NULL}, stocks)
  opts <- list(chunkSize=CHUNKSIZE, loadFactor=LOADFACTOR)
  par.time <- system.time(
    r1 <- eachElem(s, fun, rep(NUMSAMPLES, NUMTASKS), NUMSTOCKS, eo=opts))
  cat('Parallel time for', stat$numWorkers, 'workers:\n')
  print(par.time)
  cat('Shutting down the sleigh\n')
  stopSleigh(s)
}
