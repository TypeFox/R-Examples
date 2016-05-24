################################################################################
##
## $Id: tradelist.calcRanks.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests "calcRanks" method of "tradelist" class
##
################################################################################

library(portfolio)

load("tradelist.calcRanks.test.RData")

## save(tl, truth.ranks, truth.rank.sorts, truth.weights,  file = "tradelist.calcRanks.test.RData", compress = TRUE)

## tests overall behavior of calcRanks method; candidates include
## buys, sells, shorts, and covers.

tl <- portfolio:::calcRanks(tl)

## truth.ranks <- tl@ranks
## truth.rank.sorts <- tl@rank.sorts

stopifnot(
          all.equal(tl@ranks, truth.ranks),
          all.equal(tl@rank.sorts, truth.rank.sorts)
          )

## tests weighting scheme and trade-cost adjustment functions by
## assigning different weights to the different sorts.

## tl@sorts <- list(sort.input.1 = 1, sort.input.2 = 1/2)

tl.1 <- tl
tl.1@sorts <- list(sort.input.1 = 1, sort.input.2 = 1/2)

tl.1 <- portfolio:::calcRanks(tl.1)

## tests that the rank ordering produced by calcRanks is correct

stopifnot(all.equal(tl.1@ranks[order(tl.1@ranks$rank.t, decreasing = TRUE),]$id,
                    truth.weights))

## tests for missing values in sorts

## tl.2 <- portfolio:::calcRanks(tl.2)
## tl.2@ranks <- tl.2@ranks[order(row.names(tl.2@ranks)),]

## stopifnot(all.equal(row.names(tl.2@ranks), truth.row.names.2))
