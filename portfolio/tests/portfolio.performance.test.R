################################################################################
##
## $Id: portfolio.performance.test.R 384 2007-01-09 03:27:22Z enos $
##
## Tests the performance method of class portfolio.
##
################################################################################

library(portfolio)

load("portfolio.performance.test.RData")

## save(test.portfolio, empty.portfolio, test.market.data, truth, empty.truth, file = "portfolio.performance.test.RData", compress = TRUE)

result <- performance(test.portfolio, test.market.data)

empty.result <- performance(empty.portfolio, test.market.data)

## Set row.names explicitly to a character vector, since the default
## is different in R 2.5.0 than in 2.4.x and we want the test to pass
## in both.

row.names(empty.result@ret.detail) <- as.character(row.names(empty.result@ret.detail))

stopifnot(
          all.equal(result@ret, truth@ret),
          all.equal(result@profit, truth@profit),
          all.equal(result@missing.price, truth@missing.price),
          all.equal(result@missing.return, truth@missing.return),
          all.equal(result@ret.detail, truth@ret.detail),
          all.equal(result@t.plus.one, truth@t.plus.one),
          
          all.equal(empty.result, empty.truth)
          )
