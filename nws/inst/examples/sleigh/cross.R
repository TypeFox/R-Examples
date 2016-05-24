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
NUMROWS <- as.integer(Sys.getenv('NUMROWS', 500))
NUMCOLS <- as.integer(Sys.getenv('NUMCOLS', 100))
########################################################################

if (is.na(SEQUENTIAL)) {
  stop('bad value for SEQUENTIAL: ', Sys.getenv('SEQUENTIAL'))
}

x <- matrix(rnorm(NUMROWS * NUMCOLS), NUMROWS, NUMCOLS)
beta <- c(rnorm(NUMCOLS / 2, 0, 5), rnorm(NUMCOLS / 2, 0, 0.25))
y <- x %*% beta + rnorm(NUMROWS, 0, 20)
dat <- data.frame(y=y, x=x)
fold <- rep(1:10, length=NUMROWS)
fold <- sample(fold)

worker <- function(foldnumber, p) {
  fun <- function(i) {
    glmfit <- glm(y ~ ., data=dat[fold != foldnumber, 1:(i+1)])
    yhat <- predict(glmfit, newdata=dat[fold == foldnumber, 1:(i+1)])
    sum((yhat - dat[fold == foldnumber, 1]) ^ 2)
    # cat('.')
    # flush.console()
  }
  mean(sapply(1:p, fun))
}

if (SEQUENTIAL) {
  cat('Running sequentially\n')
  print(system.time(srss <- sapply(1:10, worker, NUMCOLS)))
  cat('Sequential results:', srss, '\n')
}

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
  workerinit <- function(d, f) { dat <<- d; fold <<- f }
  eachWorker(s, workerinit, dat, fold)
  print(system.time(
    prss <- unlist(eachElem(s, worker, data.frame(foldnumber=1:10, p=NUMCOLS)))))
  cat('Parallel results:', prss, '\n')
  stopSleigh(s)
}
