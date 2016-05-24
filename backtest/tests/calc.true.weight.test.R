###################################################
##
##
## Tests for function "calc.true.weight"
##
##################################################

library(backtest)

load("calc.true.weight.test.RData")

## save(calc.data, calc.truth, file = "calc.true.weight.RData")


result <- backtest:::calc.true.weight(calc.data, "date", "id", 2)
true.weight.result <- result$weight


stopifnot(
          isTRUE(all.equal(true.weight.result, calc.truth))
          )
             
