################################################################################
##
## $Id: portfolio.expandData.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests the expandData method of class portfolio.
##
################################################################################

library(portfolio)

load("portfolio.expandData.test.RData")

## save(test, truth, file = "portfolio.expandData.test.RData", compress = TRUE)

result <- expandData(test)

stopifnot(
          all.equal(truth$id, result@data$id),
          all.equal(truth$price, result@data$price),
          all.equal(truth$in.var, result@data$in.var)
          )
