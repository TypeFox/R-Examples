################################################################################
##
## $Id: calc.turnover.test.R 1300 2008-08-27 21:01:11Z zhao $
##
## Tests for function "calc.turnover"
##
################################################################################

library(backtest)

load("calc.turnover.test.RData")

## save(x.id, x.bucket, x.date, x.truth, file = "calc.turnover.test.RData", compress = TRUE)

x.result <- backtest:::calc.turnover(x.id, x.bucket, x.date)

stopifnot(
          isTRUE(all.equal(x.result, x.truth))
          )
