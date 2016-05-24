## summary.test.R
## Yang Lu Yang.Lu@williams.edu

library(pa)

## data(jan)
## b1 <- brinson(x = jan)
## truth <- summary(b1)
## data(quarter)
## b2 <- brinson(x = quarter)
## truth.multi <- summary(b2)
## save(truth, truth.multi, file = "summary.test.RData")


load("summary.test.RData")

## Single-period

data(jan)
b1 <- brinson(x = jan)
result <- summary(b1)
stopifnot(all.equal(result, truth))

## Multi-period

data(quarter)
b2 <- brinson(x = quarter)
result.multi <- summary(b2)
stopifnot(all.equal(result.multi, truth.multi))

