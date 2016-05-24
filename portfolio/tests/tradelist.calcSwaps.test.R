################################################################################
##
## $Id: tradelist.calcSwaps.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests "calcSwaps" method of "tradelist" class
##
################################################################################

library(portfolio)

load("tradelist.calcSwaps.test.RData")

## save(tl, tl.1, tl.2, tl.3, truth.swaps, truth.row.names, truth.row.names.2, truth.row.names.3, file = "tradelist.calcSwaps.test.RData", compress = TRUE)

tl  <- portfolio:::calcSwaps(tl)

stopifnot(all.equal(tl@swaps, truth.swaps))

## tests that dummy sells are created when the portfolio is
## under-invested

tl.1 <- portfolio:::calcSwaps(tl.1)

## compares generated row names to expected row names

stopifnot(all.equal(row.names(tl.1@swaps), truth.row.names))

## tests normal creation of swaps (matching buys with sells)

tl.2 <- portfolio:::calcSwaps(tl.2)

stopifnot(all.equal(row.names(tl.2@swaps), truth.row.names.2))

## tests that dummy buys are created when the portfolio is
## over-invested

tl.3 <- portfolio:::calcSwaps(tl.3)

stopifnot(all.equal(row.names(tl.3@swaps), truth.row.names.3))
