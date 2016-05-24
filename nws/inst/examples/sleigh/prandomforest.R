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
  rf <- randomForest(Class ~ ., data=SYNTH, ntree=ntree, ...)
  resp <- grep('^Class$', names(SYNTH))
  nwsStore(SleighUserNws, 'votes',
    predict(rf, SYNTH[,-resp], type='vote', norm.votes=FALSE))
  NULL
}

prandomForest <- function(s, ntree=500, ..., tasksize=ceiling(ntree/workerCount(s))) {
  numworkers <- workerCount(s)
  if (numworkers < 1) stop('no workers in sleigh')
  vntree <- rep(tasksize, floor(ntree / tasksize))
  mod <- ntree - sum(vntree)
  if (mod > 0)
    vntree <- c(vntree, mod)
  stopifnot(sum(vntree) == ntree)
  eo <- list(blocking=FALSE)
  pending <- eachElem(s, prandomForestWorker,
    list(vntree),
    list(...),
    eo=eo)
  r <- 0
  for (i in seq(along.with=vntree))
    r <- r + nwsFetch(s@userNws, 'votes')
  waitSleigh(pending)
  resp <- grep('^Class$', names(SYNTH))
  levs <- levels(SYNTH[,resp])
  factor(levs[apply(r, 1, which.max)])
}

srandomForest <- function(...) {
  rf <- randomForest(Class ~ ., data=SYNTH, ...)
  resp <- grep('^Class$', names(SYNTH))
  predict(rf, SYNTH[,-resp])
}

nalph <- max(min(NALPH, length(LETTERS)), 2)
ALPHA <- LETTERS[1:nalph]
SYNTH <- synth(NROWS, NCOLS, ALPHA)

# Sequential version
if (SEQUENTIAL) {
  cat('Running sequential random forest: ntree =', NTREE, '\n')
  t1 <- system.time(seq.pred <- srandomForest(ntree=NTREE, importance=TRUE))
  cat('Sequential time:\n')
  print(t1)
}

# Parallel version
for (numWorkers in NUMWORKERS) {
  cat('Creating sleigh\n')
  s <- if (length(NODELIST) > 0) {
         sleigh(nodeList=rep(NODELIST, length.out=numWorkers),
                launch=sshcmd, verbose=VERBOSE)
       } else {
         sleigh(workerCount=numWorkers, verbose=VERBOSE)
       }
  cat('Waiting for the sleigh to start\n')
  stat <- status(s, TRUE, 20)
  cat('Initializing parallel random forest\n')
  eachWorker(s, prandomForestInit, initfun=synth, nrows=NROWS,
             ncols=NCOLS, alpha=ALPHA)
  for (tasksize in TASKSIZE) {
    if (tasksize <= 0) tasksize <- ceiling(NTREE / numWorkers)
    cat('Running parallel random forest: numworkers =', NUMWORKERS,
        'ntree =', NTREE, 'tasksize =', tasksize, '\n')
    t2 <- system.time(
      par.pred <- prandomForest(s, ntree=NTREE, importance=TRUE, tasksize=tasksize))
    cat('Parallel time:\n')
    print(t2)

    # Check results if possible
    if (SEQUENTIAL) {
      print(all(seq.pred == par.pred))
      print(identical(seq.pred, par.pred))
    }
  }
  stopSleigh(s)
}
