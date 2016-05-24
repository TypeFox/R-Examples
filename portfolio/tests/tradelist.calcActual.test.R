################################################################################
##
## $Id: tradelist.calcActual.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests "calcActual" method of "tradelist" class
##
################################################################################

library(portfolio)

load("tradelist.calcActual.test.RData")

## save(tl, truth.actual, file = "tradelist.calcActual.test.RData", compress = TRUE)

tl <- portfolio:::calcActual(tl)

stopifnot(all.equal(tl@actual, truth.actual))

## truth.actual <- tl@actual


