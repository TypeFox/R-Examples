################################################################################
##
## $Id: portfolioSim.class.test.R 1323 2009-05-04 13:19:14Z enos $
##
## Tests the validity function for the portfolioSim class
##
################################################################################

library(portfolioSim)

## save(test.periods, file = "portfolioSim.class.test.RData", compress = TRUE)

load("portfolioSim.class.test.RData")

result.0 <- new("portfolioSim")

stopifnot(validObject(result.0))

result.1 <- new("portfolioSim",
                periods = test.periods,
                freq = 12,
                data.interface = new("sdiDf"),
                trades.interface = new("stiFromSignal"),
                start.holdings = new("portfolio"),
                fill.volume.pct = 100,
                exp.var = "blah",
                out.loc = "out/")

stopifnot(validObject(result.1))

test.periods.0 <- test.periods
names(test.periods.0)[2] <- "wrong.name"

trial.0 <- try(
               new("portfolioSim",
                   periods = test.periods.0
                   ),
               silent = TRUE
               )

stopifnot(
          all.equal(class(trial.0), "try-error"),
          as.logical(grep("Columns required.*period.*start.*end", trial.0[1]))
          )

test.periods.1 <- test.periods
test.periods.1$period <- as.complex(test.periods.1$period)

trial.1 <- try(
               new("portfolioSim",
                   periods = test.periods.1
                   ),
               silent = TRUE
               )

stopifnot(
          all.equal(class(trial.1), "try-error"),
          length(grep("column.*periods slot.*not orderable", trial.1[1])) > 0
          )

trial.2 <- try(
               new("portfolioSim",
                   periods = test.periods.1,
                   freq = c(12,1)
                   ),
               silent = TRUE
               )

stopifnot(
          all.equal(class(trial.2), "try-error"),
          as.logical(grep("freq.*numeric vector.*1", trial.2[1]))
          )
