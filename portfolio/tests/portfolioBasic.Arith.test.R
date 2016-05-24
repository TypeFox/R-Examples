################################################################################
##
## $Id: portfolioBasic.Arith.test.R 406 2007-04-19 16:30:22Z enos $
##
## Tests "+" method of "portfolioBasic"
##
################################################################################

library(portfolio)

load("portfolioBasic.Arith.test.RData")

## save(p.0, p.1, truth, file = "portfolioBasic.Arith.test.RData", compress = TRUE)

## Constructs data and portfolioBasics that should add ("+") correctly

data.0 <- data.frame(id = as.character(letters[1:20]), in.var = 1:20)
data.1 <- data.frame(id = as.character(letters[26:7]), in.var = 20:1)

p.0 <- new("portfolioBasic", id.var = "id", symbol.var = "symbol.var",
              in.var = "in.var", ret.var = "ret.var", type = "sigmoid",
              size = 10, data = data.0, sides = c("long", "short"))

p.1 <- new("portfolioBasic", id.var = "id", symbol.var = "symbol.var",
              in.var = "in.var", ret.var = "ret.var", type = "linear",
              size = 8, data = data.1, sides = c("long", "short"))

p.0@data$id <- as.character(p.0@data$id)
p.1@data$id <- as.character(p.1@data$id)

p.0@weights$id <- as.character(p.0@weights$id)
p.1@weights$id <- as.character(p.1@weights$id)

p.sum <- p.0 + p.1

stopifnot(
          isTRUE(all.equal(p.sum, truth)),
          nrow(p.sum@data) == length(union(p.0@data$id, p.1@data$id))
          )


p.pos <- p.0

p.neg <- p.0
p.neg@weights$weight <- -1 * p.neg@weights$weight

p.sum <- p.pos + p.neg

stopifnot(validObject(p.sum, test = TRUE),
          nrow(p.sum@weights) == 0)
