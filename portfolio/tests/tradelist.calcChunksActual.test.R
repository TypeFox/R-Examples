################################################################################
##
## $Id: tradelist.calcChunksActual.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests "calcChunksActual" method of "tradelist" class
##
################################################################################

library(portfolio)

load("tradelist.calcChunksActual.test.RData")

## save(tl, tl.1,  truth.chunks.actual, truth.row.names, file = "tradelist.calcChunksActual.test.RData", compress = TRUE)

tl <- portfolio:::calcChunksActual(tl)

tl.1 <- portfolio:::calcChunksActual(tl.1)
tl.1@chunks.actual <- tl.1@chunks.actual[order(row.names(tl.1@chunks.actual)),]

stopifnot(all.equal(tl@chunks.actual, truth.chunks.actual))
stopifnot(all.equal(row.names(tl.1@chunks.actual), truth.row.names))

## truth.chunks.actual <- tl@chunks.actual


