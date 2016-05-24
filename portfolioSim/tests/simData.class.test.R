################################################################################
##
## $Id: simData.class.test.R 344 2006-10-01 05:06:05Z enos $
##
## Tests the validity function for the simData class
##
################################################################################

library(portfolioSim)

## save(test.data, file = "simData.class.test.RData", compress = TRUE)

load("simData.class.test.RData")

result.0 <- new("simData", data = NULL)

result.1 <- new("simData", data = test.data[1:5,])

data.1 <- test.data[1:5,]
data.1$id <- as.factor(data.1$id)

data.2 <- test.data[1:5,]
data.2$start.price <- as.character(data.2$start.price)

data.3 <- test.data[1:5,]
data.3$end.price <- as.character(data.3$end.price)

data.4 <- test.data[1:5,]
data.4$ret <- as.character(data.4$ret)

trial.0 <- try(
               new("simData", data = test.data[1:5,1:4]), silent = TRUE
               )

trial.1 <- try(
               new("simData", data = data.1), silent = TRUE
               )

trial.2 <- try(
               new("simData", data = data.2), silent = TRUE
               )

trial.3 <- try(
               new("simData", data = data.3), silent = TRUE
               )

trial.4 <- try(
               new("simData", data = data.4), silent = TRUE
               )

trial.5 <- try(
               new("simData", data = test.data[c(-6,-8),]), silent = TRUE
               )

trial.6 <- try(
               new("simData", data = test.data[c(-6,-7),]), silent = TRUE
               )

trial.7 <- try(
               new("simData", data = test.data[c(-7,-8),]), silent = TRUE
               )

stopifnot(
          validObject(result.0),
          validObject(result.1),
          class(trial.0) == "try-error",
          class(trial.1) == "try-error",
          class(trial.2) == "try-error",
          class(trial.3) == "try-error",
          class(trial.4) == "try-error",
          class(trial.5) == "try-error",
          class(trial.6) == "try-error",
          class(trial.7) == "try-error"
          )
