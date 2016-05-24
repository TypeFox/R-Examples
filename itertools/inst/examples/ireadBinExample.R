# This is a terrible example and desperately needs to be replaced
# with something that appears somewhat useful.

library(itertools)
library(foreach)

n <- 1000
zz <- file("testbin", "wb")
expected <- foreach(1:1000) %do% {
  x <- rnorm(n)
  writeBin(x, zz)
  mean(x)
}
close(zz)

it <- ireadBin("testbin", "double", n=n)
actual <- foreach(x=it) %do% {
  mean(x)
}

print(identical(actual, expected))

unlink("testbin")
