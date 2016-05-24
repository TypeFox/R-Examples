library(foreach)
library(itertools)

n <- 10777
nrowsx <- 959
ncolsy <- 101
chunks <- 80  # applied to "n"
xchunks <- 1  # applied to "nrowsx"
ychunks <- 1  # applied to "ncolsy"

x <- matrix(rnorm(nrowsx * n), nrow=nrowsx)
y <- matrix(rnorm(n * ncolsy), nrow=n)
expected <- x %*% y

actual <-
  foreach(ia=iarray(x, c(2,1), chunks=c(chunks,xchunks)),
          .combine='rbind') %:%
    foreach(a=ia, ib=iarray(y, c(2,1), chunks=c(ychunks,chunks)),
            .combine='+') %:%
      foreach(b=ib, .combine='cbind') %do% {
        a %*% b
      }
all.equal(actual, expected)

actual <-
  foreach(ib=iarray(y, c(1,2), chunks=c(chunks,ychunks)),
          .combine='cbind') %:%
    foreach(b=ib, ia=iarray(x, c(1,2), chunks=c(xchunks,chunks)),
            .combine='+') %:%
      foreach(a=ia, .combine='rbind') %do% {
        a %*% b
      }
all.equal(actual, expected)

actual <-
  foreach(bsub=iarray(y, 2, chunks=ychunks),
          .combine='cbind') %:%
    foreach(ia=iarray(x, c(2,1), chunks=c(chunks,xchunks)),
            .combine='rbind') %:%
      foreach(a=ia, b=iarray(bsub, 1, chunks=chunks),
              .combine='+') %do% {
        a %*% b
      }
all.equal(actual, expected)

iqseq <- function(n, ...) {
  i <- 0
  it <- idiv(n, ...)

  nextEl <- function() {
    j <- i + nextElem(it)
    val <- call(':', i + 1, j)
    i <<- j
    val
  }

  obj <- list(nextElem=nextEl)
  class(obj) <- c('abstractiter', 'iter')
  obj
}

actual <-
  foreach(bcols=iqseq(ncol(y), chunks=ychunks),
          .combine='cbind') %:%
    foreach(ia=iarray(x, c(2,1), chunks=c(chunks,xchunks)),
            .combine='rbind') %:%
      foreach(a=ia, b=iarray(y, 1, chunks=chunks, idx=list(TRUE, bcols)),
              .combine='+') %do% {
        a %*% b
      }
all.equal(actual, expected)
