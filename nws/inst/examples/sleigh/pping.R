library(nws)

ping <- function (ws, totalTasks) {
  for (i in seq(length.out=totalTasks)) {
    r <- nwsFetch(ws, 'ping')
    stopifnot(length(r) == 2)
    nwsStore(ws, r[1], r[2])
  }
}

pong <- function (numTasks, size) {
  # this uses global variables SleighUserNws and SleighRank
  pong <- sprintf("pong_%d", SleighRank)
  payload <- paste(rep("#", size), collapse='')

  for (i in seq(length.out=numTasks)) {
    nwsStore(SleighUserNws, 'ping', c(pong, payload))
    j <- nwsFetch(SleighUserNws, pong)
  }
}

# read user input
cat('Please enter Host NumTasks Size (in order):\n')
input <- scan("", list(host="", numTasks=0, size=0), nmax=1, quiet=TRUE)

host <- input$host
if (length(host) == 0 || host == "") host <- "localhost"

numTasks <- as.integer(input$numTasks)
if (length(numTasks) == 0 || is.na(numTasks) || numTasks < 1)
  numTasks <- 10
  
size <- as.integer(input$size)
if (length(size) == 0 || is.na(size) || size < 1)
  size <- 10

# create a Sleigh and compute the number of workers
cat('Creating sleigh using NWS server on', host, '\n')
s <- sleigh(nwsHost=host)

numWorkers <- nwsFind(s@nws, 'workerCount')

eo <- list(blocking=FALSE)
eachWorker(s, pong, numTasks, size, eo=eo)

totalTasks <- numWorkers * numTasks
totalTime <- system.time(ping(s@userNws, totalTasks))[3]
cat("\n\nSeconds per operation: ", totalTime / (4 * totalTasks), '\n')
cat("Payload size is approximately ", size, " bytes\n\n")

stopSleigh(s)
