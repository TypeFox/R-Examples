## show.test.R
## Yang Lu Yang.Lu@williams.edu

library(pa)

## data(jan)
## b1 <- brinson(x = jan)
## truth <- b1
## save(truth, file = "show.test.RData")

load("show.test.RData")

## Single-period

data(jan)
result <- brinson(x = jan)
stopifnot(all.equal(result, truth))

