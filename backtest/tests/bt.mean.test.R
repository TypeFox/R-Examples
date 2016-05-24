################################################################################
##
## $Id: bt.mean.test.R 1300 2008-08-27 21:01:11Z zhao $
##
## Tests for function "bt.mean"
##
################################################################################

library(backtest)

load("bt.mean.test.RData")

## save(x, truth, file = "bt.mean.test.RData", compress = TRUE)

stopifnot(
          isTRUE(all.equal(backtest:::.bt.mean(x), truth))
        )
