library(nws)

ring <- function(numWorkers, numTasks) {
  mine <- sprintf("worker_%d", SleighRank)
  nextOne= (SleighRank + 1) %% numWorkers
  his <- sprintf("worker_%d", nextOne)

  for (i in 1:numTasks) {
    j <- nwsFetch(SleighUserNws, mine)
    nwsStore(SleighUserNws, his, j + 1)
  }
}

# read user input
cat("Please enter Host NumTasks (in order):\n")
input <- scan("", list(host='', numTasks=0), nmax=1, quiet=TRUE)

host <- input$host
if (length(host) == 0 || host == "")
  host <- "localhost"

numTasks <- as.integer(input$numTasks)
if (length(numTasks) == 0 || is.na(numTasks) || numTasks < 1)
  numTasks <- 10

cat('Creating sleigh using NWS server on', host, '\n')
s <- sleigh(nwsHost=host)
numWorkers <- nwsFind(s@nws, 'workerCount')

# include the master as one of the workers
numWorkers <- numWorkers + 1
cat("Number of workers (include the master):", numWorkers, "\n")

# tell the workers to execute the ring function defined in this file
eo <- list(blocking=FALSE)
eachWorker(s, ring, numWorkers, numTasks, eo=eo)

# the master becomes the last worker
SleighRank <<- numWorkers - 1
SleighUserNws <<- s@userNws

cat("Master assigned rank", SleighRank, "\n")

# time how long it takes the token to go all
# the way around the ring numTask times
totalTime <- system.time({
  nwsStore(s@userNws, 'worker_0', 0)
  ring(numWorkers, numTasks)
})[3]

token <- nwsFetch(s@userNws, 'worker_0')
stopifnot(token == numTasks * numWorkers)

cat("The token was passed", token, "times in", totalTime, "seconds\n")
cat("Seconds per operation:", totalTime/ (2 * token), "\n")

stopSleigh(s)
