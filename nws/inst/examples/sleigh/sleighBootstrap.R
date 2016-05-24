# parallel bootstrapping

init <- function() {
  nwsDir <- .path.package('nws', quiet=TRUE)
  source(file.path(nwsDir, 'examples', 'sleigh', 'nuclearBootstrapInit.R'))
}

task <- function(y)
  boot(nuke.data, nuke.fun, R=y, m=1, fit.pred=new.fit, x.pred=new.data)

library(nws)
s <- sleigh()

R <- 20000
numChunks <- 200
chunkSize <- R / numChunks
t <- eachWorker(s, init)
stime <- system.time(eachElem(s, task, rep(chunkSize, numChunks)))
cat('Elapsed time:', stime[3], '\n')
