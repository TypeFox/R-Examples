################################################################################
##
## $Id: portfolio.mvLong.test.R 410 2007-04-22 00:13:06Z enos $
##
## Tests for method "mvLong".
##
################################################################################

library(portfolio)

price <- rnorm(20, mean = 80, sd = 30)

data.0 <- data.frame(id = 1:20, symbol.var = letters[1:20],
                     in.var = price,
                     ret.var = rnorm(20, mean = 0, sd = 0.1),
                     price.var = price, round.lot = 1)

## corner case of 0 positions

p <- new("portfolio", id.var = "id", symbol.var = "symbol.var",
         in.var = "in.var", ret.var = "ret.var", type = "equal",
         size = "quintile", equity = 100000, sides = "long",
         price.var = "price.var", data = data.0)

p@shares <- p@shares[0,]

stopifnot(portfolio:::mvLong(p) == 0)
