## regress.test.R
## Yang Lu Yang.Lu@williams.edu

library(pa)

## data(jan)
## truth <- regress(x = jan)
## save(truth, file = "regress.test.RData")

load("regress.test.RData")

## Single-period

data(jan)
result <- regress(x = jan)
stopifnot(all.equal(result, truth))


