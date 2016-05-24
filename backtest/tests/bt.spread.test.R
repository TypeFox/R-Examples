################################################################################
##
## $Id: bt.spread.test.R 1300 2008-08-27 21:01:11Z zhao $
##
## Tests for function "bt.spread"
##
################################################################################

library(backtest)

load("bt.spread.test.RData")

## save(m, n, sd, truth, file = "bt.spread.test.RData", compress = TRUE)

stopifnot(
          all(mapply(all.equal, backtest:::.bt.spread(m, n, sd), truth))
          )
