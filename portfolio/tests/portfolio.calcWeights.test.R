################################################################################
##
## $Id: portfolio.calcWeights.test.R 1625 2010-02-18 19:44:29Z enos $
##
## Tests "calcWeights" method of "portfolio"
##
################################################################################

library(portfolio)

## Invalid Input: "price.var" must be in names(x)

## constructs a usable data.frame and portfolio

data.0 <- data.frame(id = 1:20, symbol.var = letters[1:20], in.var = 1:20,
                     ret.var = 1:20, price.var = 1:20)
p.0 <- new("portfolio", id.var = "id", symbol.var = "symbol.var",
                   in.var = "in.var", ret.var = "ret.var", type = "equal",
                   size = "quintile", equity = 100000,
                   price.var = "price.var", data = data.0)

p.0@price.var <- "foo.var"

trial.0 <- try(calcWeights(p.0), silent = TRUE)

if(class(trial.0) == "try-error"){
  stopifnot(as.logical(grep("valid price.var!", trial.0[1])))
}

p.0@price.var <- "price.var"

## Corner Case: shares slot has 0 rows

p.0@shares <- p.0@shares[0,]

## "calcWeights" returns "object" under this condition

stopifnot(class(calcWeights(p.0)) == "portfolio")

## Tests misuse of id.var and the "id" column when
## id.var != "id"

## constructs a usable data.frame and portfolio

data.1 <- data.frame(ric = 1:20, symbol.var = letters[1:20], in.var = 1:20,
                     ret.var = 1:20, price.var = 1:20)

p.2 <- new("portfolio", id.var = "ric", symbol.var = "symbol.var",
              in.var = "in.var", ret.var = "ret.var", type = "equal",
              size = "quintile", equity = 100000,
              price.var = "price.var", data = data.1)

## makes portfolio invalid by changing the name of "id" column

names(p.2@data)[pmatch("id", names(p.2@data))] <- "ric"
trial.1 <- try(calcShares(p.2), silent = TRUE)

if(class(trial.1) == "try-error"){
  stopifnot(isTRUE(as.logical(grep("\\'by\\'.*valid.column", trial.1[1]))))
}
