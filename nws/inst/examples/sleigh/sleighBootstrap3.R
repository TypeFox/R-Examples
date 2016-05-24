# parallel bootstrapping

source('nuclearBootstrapInit.R')

library(nws)
s <- sleigh()

R <- 20000
y <- ceiling(R / workerCount(s))  # rounding up
f <- function(...) { library(boot); boot(...) }
results <- eachWorker(s, f, nuke.data, nuke.fun, R=y, m=1, fit.pred=new.fit, x.pred=new.data)

for (r in results) {
  get(getOption("device"))()
  plot(r)
}
