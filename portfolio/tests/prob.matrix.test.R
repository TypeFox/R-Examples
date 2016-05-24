################################################################################
##
## $Id: prob.matrix.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests that ".prob.matrix" correctly forms probabilities from a matrix of
## distances
##
################################################################################

library(portfolio)

## save(truth, file = "prob.matrix.test.RData", compress = TRUE)

load("prob.matrix.test.RData")

## tests that the match.matrix matches the values with the highest
## probabilities

## sets all the values of i to a low probability

set.seed(1)
distances <- matrix(abs(rnorm(30, sd = .30)),
                    nrow = 5,
                    ncol = 6,
                    dimnames = list(1:5, letters[1:6])
            )

## sets letter "a" to have a high probability of being chosen

test <- portfolio:::.prob.matrix(distances)

stopifnot(
          all.equal(truth, test)
          )
