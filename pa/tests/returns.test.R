## returns.test.R
## Yang Lu Yang.Lu@williams.edu

library(pa)

## data(jan)
## b1 <- brinson(x = jan)
## truth <- returns(b1)

## r1 <- regress(jan)
## truth.r1 <- returns(r1)

## data(quarter)
## b2 <- brinson(x = quarter)
## truth.multi <- returns(b2)

## r2 <- regress(quarter)
## truth.multi.r2 <- returns(r2)

## save(truth, truth.multi, truth.r1, truth.multi.r2, file = "returns.test.RData")

load("returns.test.RData")

## Single-period
data(jan)
b1 <- brinson(x = jan)
result <- returns(b1, var = "sector")
stopifnot(all.equal(result, truth))

r1 <- regress(jan)
result.r1 <- returns(r1, var = "sector")
stopifnot(all.equal(result.r1, truth.r1))

## Multi-period

data(quarter)
b2 <- brinson(x = quarter)
result.multi <- returns(b2, var = "sector")
stopifnot(all.equal(result.multi, truth.multi))

r2 <- regress(quarter)
result.multi.r2 <- returns(r2, var = "sector")
stopifnot(all.equal(result.multi.r2, truth.multi.r2))

