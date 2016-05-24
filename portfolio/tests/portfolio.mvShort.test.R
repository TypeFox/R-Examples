################################################################################
##
## $Id: portfolio.mvShort.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests for method "mvShort".
##
################################################################################

library(portfolio)

price <- rnorm(20, mean = 80, sd = 30)

data.1 <- data.frame(id = 1:20, symbol.var = letters[1:20],
                     in.var = price,
                     ret.var = rnorm(20, mean = 0, sd = 0.1),
                     price.var = price, round.lot = 1)

## corner case of no positions

portfolio.1 <- new("portfolio", id.var = "id", symbol.var = "symbol.var",
               in.var = "in.var", ret.var = "ret.var", type = "equal",
               size = "quintile", equity = 100000, sides = "short",
               price.var = "price.var", data = data.1)

portfolio.1@shares <- portfolio.1@shares[0,]

stopifnot(portfolio:::mvShort(portfolio.1) == 0)
