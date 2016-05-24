##############################################################
##
## Another test for function "bucketize"
##
##############################################################

library(backtest)

load("bucketize.test2.RData")


## save(tmp.2, tmp.2.x, tmp.2.y, truth.2, file = "bucketize.test2.RData")

## Find result from a bucketize test

result.2 <- backtest:::bucketize(tmp.2, tmp.2.x, tmp.2.y, compute = mean)

## Did the test case work?

stopifnot(
          isTRUE(all.equal(result.2, truth.2))
          )
