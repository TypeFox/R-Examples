library(nws)
library(randomForest)

as.charlist <- function(s) {
  x <- strsplit(s, ' +', extended=TRUE)[[1]]
  x[x != '']
}

as.intlist <- function(s) {
  x <- as.integer(strsplit(s, ' +', extended=TRUE)[[1]])
  x[!is.na(x)]
}

######################## Configuration Section #########################
NODELIST <- as.charlist(Sys.getenv('NODELIST', ''))
NUMWORKERS <- if (length(NODELIST) > 0) length(NODELIST) else 4
NUMWORKERS <- as.intlist(Sys.getenv('NUMWORKERS', NUMWORKERS))
SEQUENTIAL <- as.logical(Sys.getenv('SEQUENTIAL', 'FALSE'))
VERBOSE <- as.logical(Sys.getenv('VERBOSE', 'FALSE'))
NROWS <- as.integer(Sys.getenv('NROWS', 100))
NCOLS <- as.integer(Sys.getenv('NCOLS', 4))
NALPH <- as.integer(Sys.getenv('NALPH', 16))
NTREE <- as.integer(Sys.getenv('NTREE', 500))
TASKSIZE <- as.intlist(Sys.getenv('TASKSIZE', 0))
########################################################################

if (is.na(SEQUENTIAL)) {
  stop('bad value for SEQUENTIAL: ', Sys.getenv('SEQUENTIAL'))
}

synth <- function(nrows, ncols, alpha=LETTERS, seed=7442) {
  if (seed > 0) set.seed(seed)
  ncols <- min(ncols, length(alpha))
  fun <- function(i, ...) sample(...)
  cols <- lapply(seq(length.out=ncols + 1), fun, alpha, nrows, replace=TRUE)
  df <- do.call('data.frame', cols)
  names(df) <- c(alpha[1:ncols], 'Class')
  df
}

prandomForestInit <- function(initfun, ...) {
  library(randomForest)
  SYNTH <<- initfun(...)
  NULL
}

prandomForestWorker <- function(ntree, ...) {
  randomForest(Class ~ ., data=SYNTH, ntree=ntree, ...)
}

prandomForest <- function(s, ntree=500, ..., tasksize=ceiling(ntree/workerCount(s))) {
  numworkers <- workerCount(s)
  if (numworkers < 1) stop('no workers in sleigh')
  vntree <- rep(tasksize, floor(ntree / tasksize))
  mod <- ntree - sum(vntree)
  if (mod > 0)
    vntree <- c(vntree, mod)
  stopifnot(sum(vntree) == ntree)
  eo <- list(blocking=TRUE)
  results <- eachElem(s, prandomForestWorker,
    list(vntree),
    list(...),
    eo=eo)
  old.warn <- options(warn=-1)$warn
  on.exit(options(warn=old.warn))
  do.call('combine', results)
}

srandomForest <- function(...) {
  randomForest(Class ~ ., data=SYNTH, ...)
}

nalph <- max(min(NALPH, length(LETTERS)), 2)
ALPHA <- LETTERS[1:nalph]
SYNTH <- synth(NROWS, NCOLS, ALPHA)

# Sequential random forest
if (SEQUENTIAL) {
  cat('Running sequential random forest: ntree =', NTREE, '\n')
  t1 <- system.time(seq.rf <- srandomForest(ntree=NTREE,
    importance=TRUE, norm.votes=FALSE))
  cat('Sequential time:\n')
  print(t1)
}

# Parallel random forest
for (numWorkers in NUMWORKERS) {
  cat('Creating a sleigh\n')
  s <- if (length(NODELIST) > 0) {
         sleigh(nodeList=rep(NODELIST, length.out=numWorkers),
                launch=sshcmd, verbose=VERBOSE)
       } else {
         sleigh(workerCount=numWorkers, verbose=VERBOSE)
       }
  cat('Waiting for the sleigh to start\n')
  stat <- status(s, TRUE, 60)
  cat('Initializing parallel random forest\n')
  eachWorker(s, prandomForestInit, initfun=synth, nrows=NROWS,
             ncols=NCOLS, alpha=ALPHA)
  for (tasksize in TASKSIZE) {
    if (tasksize <= 0) tasksize <- ceiling(NTREE / numWorkers)
    cat('Running parallel random forest: numworkers =', numWorkers,
        'ntree =', NTREE, 'tasksize =', tasksize, '\n')
    t2 <- system.time(
      par.rf <- prandomForest(s, ntree=NTREE, importance=TRUE,
        norm.votes=FALSE, tasksize=tasksize))
    cat('Parallel time:\n')
    print(t2)
  }
  stopSleigh(s)
}
