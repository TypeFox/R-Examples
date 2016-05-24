## brinson.test.R
## Yang Lu Yang.Lu@williams.edu

library(pa)

## data(jan)
## truth <- brinson(x = jan)
## save(truth, file = "brinson.test.RData")

load("brinson.test.RData")

## Single-period

data(jan)
result <- brinson(x = jan)

stopifnot(all.equal(result, truth))

