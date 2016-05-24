# parallel bootstrapping

task <- function(y) {
  nwsDir <- .path.package('nws', quiet=TRUE)
  source(file.path(nwsDir, 'examples', 'sleigh', 'nuclearBootstrapInit.R'))
  boot(nuke.data, nuke.fun, R=y, m=1, fit.pred=new.fit, x.pred=new.data)
}

library(nws)
s <- sleigh()

R <- 20000
chunkSize <- ceiling(R / workerCount(s))  # rounding up
results <- eachWorker(s, task, chunkSize)

library(boot)
for (r in results) {
  get(getOption("device"))()
  plot(r)
}
