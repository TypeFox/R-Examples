# simple program that multiplies value by 2.
library(nws)

worker <- function() {
  repeat {
    task <- nwsFetch(SleighUserNws, 'task')
    result <- 2 * task
    nwsStore(SleighUserNws, 'result', c(task, result))
  }
}

# change launch if you add nodeList parameter
s <- sleigh()
eo <- list(blocking=FALSE)
eachWorker(s, worker, eo=eo)

# read user input
input <- readline('Please enter number of tasks:\n')
numTasks <- as.integer(input)
if (is.na(numTasks) || numTasks < 1)
  stop("Please enter a positive number")

for (i in seq(length.out=numTasks))
  nwsStore(s@userNws, 'task', i)

for (i in seq(length.out=numTasks)) {
  result <- nwsFetch(s@userNws, 'result')
  cat(result[1], 'times 2 is', result[2], '\n')
}
  
stopSleigh(s)
