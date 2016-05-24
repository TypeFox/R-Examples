# parallel bootstrapping

pboot <- function(s, R=20000) {
  source('nuclearBootstrapInit.R', local=TRUE)
  chunkSize <- ceiling(R / workerCount(s))  # rounding up
  task <- function(y) {
    library(boot)
    boot(nuke.data, nuke.fun, R=y, m=1, fit.pred=new.fit, x.pred=new.data)
  }
  eachWorker(s, task, chunkSize)
}

library(nws)
s <- sleigh()

results <- pboot(s)

for (r in results) {
  get(getOption("device"))()
  plot(r)
}
