###################################################
##
##
## Tests for function "overlaps.compute"
##
##################################################

library(backtest)

load("overlaps.compute.test.RData")

## save(over.data, weight.truth, file = "overlaps.compute.test.RData")


result <- backtest:::overlaps.compute(over.data, "in.factor", "date", "id", 2)
weight.result <- result$weight


stopifnot(
          isTRUE(all.equal(weight.truth, weight.result))
          )


