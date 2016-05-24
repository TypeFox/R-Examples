################################################################################
##
## $Id: portfolio.Arith.test.R 406 2007-04-19 16:30:22Z enos $
##
## Tests for adding portolios
##
################################################################################

library(portfolio)

load("portfolio.Arith.test.RData")

## save(p.0, p.1, truth, file = "portfolio.Arith.test.RData", compress = TRUE)

p.sum <- p.0 + p.1

stopifnot(
          isTRUE(all.equal(p.sum, truth)),
          nrow(p.sum@data) == length(union(p.0@data$id, p.1@data$id))
          )

p.pos <- p.0

p.neg <- p.0
p.neg@shares$shares <- -1 * p.neg@shares$shares
p.neg <- calcWeights(p.neg)

p.sum <- p.pos + p.neg

stopifnot(validObject(p.sum, test = TRUE),
          nrow(p.sum@shares) == 0)
