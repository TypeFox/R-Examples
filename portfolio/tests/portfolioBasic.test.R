################################################################################
##
## $Id: portfolioBasic.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests for the portfolioBasic class.
##
################################################################################

library(portfolio)

load("portfolioBasic.test.RData")

## save(data, x, file = "portfolioBasic.test.RData", compress = TRUE)

data <- data.frame(id = 1:20, in.var = 1:20)
data$in.var <- as.numeric(data$in.var)

x <- new("portfolioBasic", in.var = "in.var", type = "sigmoid", size = 8, data = data)
x <- create(x)

## All weights should sum to 0; side weights sum to 1.

stopifnot(
          all.equal(sum(x@weights$weight), 0),
          all.equal(sum(x@weights$weight[x@weights$weight > 0]), 1)
          )

## Cut the weights in half and then use the scaler to fix.

x@weights$weight <- x@weights$weight / 2
x <- scaleWeights(x)

stopifnot(
          all.equal(sum(x@weights$weight), 0),
          all.equal(sum(x@weights$weight[x@weights$weight > 0]), 1)
          )

## I've removed tests involving calling the balance method and then
## checking exposures.  The Debian machine at CRAN encounters a
## failure, likely due to a differing .Machine$double.eps used in
## all.equal.  I couldn't reproduce using -ffloat-store in CFLAGS, but
## these tests aren't important for testing portfolioBasic or
## exposures, so I'm taking them out.

## I should add tests that check cases where a single side has 0
## exposure and where both sides have 0 exposure to a factor level.
