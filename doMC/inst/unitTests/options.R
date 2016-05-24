test01 <- function() {
  x <- list(1:3, 1:9, 1:19)
  cs <- 1:20

  for (chunkSize in cs) {
    mcopts <- list(preschedule=FALSE)
    for (y in x) {
      actual <- foreach(i=y, .options.multicore=mcopts) %dopar% i
      checkEquals(actual, as.list(y))
      actual <- foreach(i=y, .combine='c', .options.multicore=mcopts) %do% i
      checkEquals(actual, y)
    }
  }
}
