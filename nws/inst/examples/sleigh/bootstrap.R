library(nws)
library(boot)

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
NUMREPS <- as.integer(Sys.getenv('NUMREPS', 20000))
NUMCHUNKS <- as.integer(Sys.getenv('NUMCHUNKS', 200))
########################################################################

if (is.na(SEQUENTIAL)) {
  stop('bad value for SEQUENTIAL: ', Sys.getenv('SEQUENTIAL'))
}

# bootstrap statistic function
nuke.fun <- function(dat, inds, i.pred, fit.pred, x.pred) {
  assign('.inds', inds, envir=.GlobalEnv)
  lm.b <- glm(fit+resid[.inds] ~ date+log(cap)+ne+ct+log(cum.n)+pt, data=dat)
  pred.b <- predict(lm.b, x.pred)
  remove('.inds', envir=.GlobalEnv)
  c(coef(lm.b), pred.b - (fit.pred + dat$resid[i.pred]))
}

# worker init function
initWorker <- function() {
  library(boot)
  assign('boot.data', nwsFind(SleighUserNws, 'boot.data'), envir=globalenv())
  NULL
}

# worker clean up function
cleanWorker <- function() {
  rm('boot.data', envir=globalenv())
}

# worker task function
pbootWorker <- function(repNo, repCounts,
                        statistic, sim, stype,
                        strata, L, m, weights,
                        ran.gen, mle, ...)
{
  nreps <- repCounts[repNo + 1] - repCounts[repNo]
  boot(boot.data, statistic, nreps, sim=sim, stype=stype, strata=strata,
       L=L, m=m, weights=weights, ran.gen=ran.gen, mle=mle, ...)
}

# parallel bootstrap function, almost like the standard boot function
pboot <- function(s, nchunks,
                  data, statistic, R, sim='ordinary', stype='i',
                  strata=rep(1, n), L=NULL, m=0, weights=NULL,
                  ran.gen=function(d, p) d, mle=NULL, ...)
{
  thisCall <- match.call()
  repCounts <- floor(R * 0:nchunks / nchunks)
  n <- if (length(dim(data)) == 2) nrow(data) else length(data)

  nwsStore(s@userNws, 'boot.data', data)
  eachWorker(s, initWorker)
  nwsDeleteVar(s@userNws, 'boot.data')

  res <- eachElem(s, pbootWorker, list(1:nchunks),
          list(repCounts, statistic, sim, stype, strata,
          L, m, weights, ran.gen, mle, ...))

  eachWorker(s, cleanWorker)

  out <- res[[1]]
  if (! is.null(out) && class(out) == 'boot') {
    t <- lapply(res, '[[', 't')
    out$t <- do.call('rbind', t)
    out$R <- R
    out$call <- thisCall
    class(out) <- c('pboot', 'boot')
  } else {
    warning('pboot returning object that is not of class "boot"', immediate.=TRUE)
  }
  out
}

# initialize the data needed to do the bootstrapping
nuke <- nuclear[,c(1,2,5,7,8,10,11)]
nuke.lm <- glm(log(cost) ~ date+log(cap)+ne+ct+log(cum.n)+pt, data=nuke)
nuke.diag <- glm.diag(nuke.lm)
nuke.res <- nuke.diag$res * nuke.diag$sd
nuke.res <- nuke.res - mean(nuke.res)
nuke.data <- data.frame(nuke, resid=nuke.res, fit=fitted(nuke.lm))
new.data <- data.frame(cost=1, date=73.00, cap=886, ne=0, ct=0, cum.n=11, pt=1)
new.fit <- predict(nuke.lm, new.data)
R <- NUMREPS

# sequential bootstrap
if (SEQUENTIAL) {
  cat('Running sequentially\n')
  seq.time <- system.time(
    seq.boot <- boot(nuke.data, nuke.fun, R=R, m=1, fit.pred=new.fit, x.pred=new.data))
  cat('Sequential time:\n')
  print(seq.time)
}

# parallel bootstrap
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
  cat('Running on', stat$numWorkers, 'workers\n')
  par.time <- system.time(
    par.boot <- pboot(s, NUMCHUNKS, nuke.data, nuke.fun, R=R, m=1, fit.pred=new.fit, x.pred=new.data))
  cat('Parallel time for', stat$numWorkers, 'workers:\n')
  print(par.time)
  cat('Shutting down the sleigh\n')
  stopSleigh(s)
}
