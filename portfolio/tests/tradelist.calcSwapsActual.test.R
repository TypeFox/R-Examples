################################################################################
##
## $Id: tradelist.calcSwapsActual.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests "calcSwapsActual" method of "tradelist" class
##
################################################################################

library(portfolio)

load("tradelist.calcSwapsActual.test.RData")

## save(tl, tl.1, tl.2, truth.swaps.actual, truth.row.names,  file = "tradelist.calcSwapsActual.test.RData", compress = TRUE)

## tests that overall behavior of method is correct

tl <- portfolio:::calcSwapsActual(tl)

## truth.swaps.actual <- tl@swaps.actual

stopifnot(all.equal(tl@swaps.actual, truth.swaps.actual))

## tests that setting turnover to 0 causes all swaps to be removed

tl.1 <- portfolio:::calcSwapsActual(tl.1)

stopifnot(all.equal(nrow(tl.1@swaps.actual), 0))

## tests that setting rank.gain.min to 3 removes all swaps with a
## rank.gain of less than 3

tl.2 <- portfolio:::calcSwapsActual(tl.2)

stopifnot(all.equal(row.names(tl.2@swaps.actual), truth.row.names))


