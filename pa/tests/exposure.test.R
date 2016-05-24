## exposure.test.R
## Yang Lu Yang.Lu@williams.edu

library(pa)

## data(jan)
## b1 <- brinson(x = jan)
## truth <- exposure(b1, var = "sector")
## data(quarter)
## b2 <- brinson(x = quarter)
## truth.multi <- exposure(b2, var = "sector")
## save(truth, truth.multi, file = "exposure.test.RData")

load("exposure.test.RData")

## Single-period

data(jan)
b1 <- brinson(x = jan)
result <- exposure(b1, var = "sector")
stopifnot(all.equal(result, truth))

## Multi-period

data(quarter)
b2 <- brinson(x = quarter)
result.multi <- exposure(b2, var = "sector")
stopifnot(all.equal(result.multi, truth.multi))

