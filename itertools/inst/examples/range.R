library(itertools)
library(foreach)
library(abind)

a <- array(rnorm(60), c(5,4,3))
expected <- apply(a, c(2,3), range)
acomb <- function(...) abind(..., along=3)

actual <-
  foreach(ia=iarray(a, c(2,3)), .combine='acomb', .multicombine=TRUE) %:%
    foreach(x=ia, .combine='cbind') %do%
      range(x)
dimnames(actual) <- NULL
print(identical(actual, expected))

actual <-
  foreach(x=iarray(a, 3, chunks=2), .combine='acomb', .multicombine=TRUE) %do%
    apply(x, c(2,3), range)
dimnames(actual) <- NULL
print(identical(actual, expected))
