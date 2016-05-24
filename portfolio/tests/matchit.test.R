################################################################################
##
## $Id: matchit.test.R 1313 2008-10-31 19:24:34Z enos $
##
## Tests for the matching method of portfolioBasic
##
################################################################################

library(portfolio)

## save(truth, truth.2, truth.3, truth.4, file = "matchit.test.RData", compress = TRUE)

load("matchit.test.RData")

data(assay)

x <- assay
x <- assay[assay$country == "USA", c("symbol", "name", "sector", "liq", "on.fl")]

## universe for test case includes all US stocks, 10 from the focus
## list, 10 identified as good matches by the matchit method, and 10
## other US stocks

all.stocks <- c("76143", "18027", "14730", "6961", "6930", "69571", "71262",
                "21266", "7308", "11746", "27043", "37495", "74206", "79463", "2923", "8267",
                "33105", "26322", "68150", "71570", "22101", "19167", "39252", "13776",
                "83265", "71301", "7631", "29780", "3604", "28225")

x <- x[all.stocks,]

for(i in names(x)){
  if(is.factor(x[[i]])){
    x[[i]] <- as.character(x[[i]])
  }
}

## done preparing data, tests greedy algorithm of entire "matchit"
## function.

test <- portfolio:::.matchit(on.fl ~ sector + liq, data = x)

## truth is a matrix

stopifnot(
          all.equal(dimnames(test)[1], dimnames(truth)[1]),
          all(mapply(all.equal, test, truth))
          )

## corner case: number of controls is less than the number of treated.

x.sub <- x[2:15,]

test.3 <- portfolio:::.matchit(on.fl ~ sector + liq, data = x.sub)

stopifnot(
          all.equal(dimnames(test.3)[1], dimnames(truth.3)[1]),
          all.equal(test.3[1,], truth.3[1,])
          )


################################################################################
## Exact
################################################################################

test <- portfolio:::.matchit(on.fl ~ sector, data = x, method = "greedy", exact = "sector")
test <- test[sort(row.names(test)), , drop = FALSE]

stopifnot(all.equal(test, truth.4),
          all.equal(table(x$sector[x$on.fl]), table(x[test,"sector"]))
          )

################################################################################
## Random matching, with exact
################################################################################

set.seed(1)

x$foo.1 <- c("a","a","b")
x$foo.2 <- c("a","b")

test <- portfolio:::.matchit(treat.var = "on.fl",  data = x, method = "random", n.matches = 10, exact = c("foo.1","foo.2"))

x.treat <- x[x$on.fl,]
for(i in 1:ncol(test)){
  x.matched <- x[test[,1],]
  stopifnot(all.equal(table(x.treat$foo.1, x.treat$foo.2),
                      table(x.matched$foo.1, x.matched$foo.2)))
}

################################################################################
## Matchit subfunctions
################################################################################

## tests that .valid.data converts character vectors to factors

x <- portfolio:::.data.prep(treat.var = "on.fl", c("sector", "country"), x)

stopifnot(is.factor(assay$sector),
          is.factor(assay$country))

