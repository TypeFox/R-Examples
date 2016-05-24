################################################################################
##
## $Id: portfolio.calcShares.test.R 1625 2010-02-18 19:44:29Z enos $
##
## Tests "calcShares" method of "portfolio"
##
################################################################################

library(portfolio)

load("portfolio.calcShares.test.RData")

## save(truth, file = "portfolio.calcShares.test.RData", compress = TRUE)

## constructs a usable data.frame and portfolio

data.0 <- data.frame(id = 1:20, symbol.var = letters[1:20], in.var = 1:20,
                     ret.var = 1:20, price.var = 1:20)
p.0 <- new("portfolio", id.var = "id", symbol.var = "symbol.var",
                   in.var = "in.var", ret.var = "ret.var", type = "equal",
                   size = "quintile", equity = 100000,
                   price.var = "price.var", data = data.0)

## tests output against precalculated values

shares <- calcShares(p.0)@shares
shares <- shares[order(as.numeric(shares$id)),]
stopifnot(all.equal(truth, shares$shares))

## Invalid Input: "price.var" not in names(data)

p.0@price.var <- "foo"

trial.0 <- try(
               calcShares(p.0),
               silent = TRUE
               )

if(class(trial.0) == "try-error"){
  stopifnot(isTRUE(as.logical(grep("Error.*valid price\\.var!", trial.0[1]))))
}

## Corner Case: shares slot has 0 rows

p.0@weights <- p.0@weights[0,]

## "calcShares" returns "object" when shares slot
## has 0 rows

stopifnot(class(calcShares(p.0)) == "portfolio")
