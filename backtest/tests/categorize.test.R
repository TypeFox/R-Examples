################################################################################
##
## $Id: categorize.test.R 1300 2008-08-27 21:01:11Z zhao $
##
## Tests for function "categorize"
##
################################################################################

library(backtest)

load("categorize.test.RData")

## save(tmp.1, tmp.1.n, truth.1, tmp.2, truth.2, file = "categorize.test.RData", compress = TRUE)

result.1 <- backtest:::categorize(tmp.1, n = tmp.1.n)
result.2 <- backtest:::categorize(tmp.2)

stopifnot(
          isTRUE(all.equal(result.1, truth.1)),
          isTRUE(all.equal(result.2, truth.2))
        )
