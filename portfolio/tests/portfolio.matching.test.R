################################################################################
##
## $Id: portfolio.matching.test.R 389 2007-01-10 04:28:44Z enos $
##
## 
##
################################################################################

library(portfolio)

load("portfolio.matching.test.RData")
## save(p, p.truth, file = "portfolio.matching.test.RData", compress = TRUE)

p.m <- matching(p, covariates = c("country", "sector", "liquidity"))
p.test <- portfolio:::.match.as.portfolioBasic(p.m, 1)

stopifnot(
          validObject(p.m),
          all.equal(p.test, p.truth)
          )

## basic test of "sample" method

p.m <- matching(p, covariates = c("sector", "liquidity"), method = "sample",
                n.matches = 5)


stopifnot(
          all.equal(dim(p.m@matches), c(33, 5))
          )

set.seed(1)

p.m <- matching(p, method = "random", n.matches = 5)


stopifnot(
          all.equal(dim(p.m@matches), c(33, 5))
          )

################################################################################
## Subroutine tests 
################################################################################

## .matching.prep

#test <- portfolio:::.matching.prep(data = p@data, weights = p@weights,
#                       covariates = c("sector", "liquidity"))

#stopifnot(
#          all(test[test$treatment, "id"] %in% p@weights$id)
#          )

## .matching.scaled.weights

id.map <- matrix(nrow = 31,
                 ncol = 1,
                 dimnames = list(p@weights$id[-(1:2)], 1)
                 )


test <- portfolio:::.matching.scale.weights(weights = p@weights, id.map = id.map)

stopifnot(
          all.equal(test$weight, rep(-0.032, length(test$weight)),
                    tolerance = 0.01)
          )

## tests the "calc.scaling.factor.R" function

orig.weights <- rep(c(0.2, -0.2), length.out = 10)
matched.weights <- rep(c(0.1, -0.1), length.out = 10)

test  <- portfolio:::.calc.scaling.factor(orig.weights, matched.weights)
truth <- c(2,2)
names(truth) <- c("-1", "1")

stopifnot(
          all.equal(test, truth)
          )

## tests the ".scale.weights" function

## long-only portfolio

scaling.factors <- 5
names(scaling.factors) <-  "1"

x <- rep(0.04, length.out = 5)

scaled.x <- portfolio:::.scale.weights(x, scaling.factors)

stopifnot(
          all.equal(sum(scaled.x), 1)
          )

## short-only portfolio

scaling.factors <- 5
names(scaling.factors) <-  "-1"

x <- rep(-0.04, length.out = 5)

scaled.x <- portfolio:::.scale.weights(x, scaling.factors)

stopifnot(
          all.equal(sum(scaled.x), -1)
          )

## long-short portfolio

scaling.factors <- c(5,5)
names(scaling.factors) <- c("-1", "1")

x <- rep(c(-0.04, 0.04), length.out = 10)

scaled.x <- portfolio:::.scale.weights(x, scaling.factors)

stopifnot(
          all.equal(sum(scaled.x), 0)
          )

## corner case where original portfolio is long-short and matched
## portfolio is long or short only

scaling.factors <- c(2,2)
names(scaling.factors) <- c("-1", "1")

x <- rep(0.1, length.out = 5)

scaled.x <- portfolio:::.scale.weights(x, scaling.factors)

stopifnot(
          all.equal(sum(scaled.x), 1)
          )
