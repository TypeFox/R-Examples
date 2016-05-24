################################################################################
##
## $Id: sdiDf.getSimData.test.R 372 2006-10-04 13:30:42Z enos $
##
## Tests for the getSimData method of class sdiDf
##
################################################################################

library(portfolioSim)

load("sdiDf.getSimData.test.RData")

## save(test.1, test.2, test.3, test.4, truth.1, truth.2, truth.3, truth.4, file = "sdiDf.getSimData.test.RData", compress = TRUE)

result.1 <- getSimData(test.1, period = "Y", verbose = FALSE)
result.2 <- getSimData(test.2, period = as.Date("1995-01-01"), verbose = FALSE)
result.3 <- getSimData(test.3, period = as.POSIXct("1994-12-31 19:00", tz = "GMT"), verbose = FALSE)
result.4 <- getSimData(test.4, period = FALSE, verbose = FALSE)

stopifnot(
          all.equal(result.1, truth.1),
          all.equal(result.2, truth.2),
          all.equal(result.3, truth.3),
          all.equal(result.4, truth.4)
          )
