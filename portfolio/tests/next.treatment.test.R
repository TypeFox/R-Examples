################################################################################
##
## $Id: next.treatment.test.R 1313 2008-10-31 19:24:34Z enos $
##
## Tests for the matching method of portfolioBasic
##
################################################################################

library(portfolio)

## save(truth, file = "next.treatment.test.RData", compress = TRUE)

load("next.treatment.test.RData")

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
x$matches <- NA
x$ps <- fitted(glm(on.fl ~ sector + liq, x, family = binomial("logit")))

treatments <- character()

while(any(is.na(x[x[["on.fl"]] == TRUE, "matches"]))){
  
  current.treatment <- portfolio:::.next.treatment(x, "on.fl")
  x[current.treatment, "matches"] <- 1
  treatments <- append(treatments, row.names(x[current.treatment,]))
  treatments
    }

stopifnot(
          isTRUE(all.equal(treatments, truth))
          )
