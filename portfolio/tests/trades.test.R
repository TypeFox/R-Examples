################################################################################
##
## $Id: trades.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests the validity function for class "trades"
##
################################################################################

library(portfolio)

load("trades.test.RData")

## save(test.data, file = "trades.test.RData", compress = TRUE)

t.base <- new("trades", trades = test.data)

stopifnot(validObject(t.base))

## Multiple ids for same side

t.1 <- t.base
t.1@trades[5,2] <- "C"

trial.1 <- try(validObject(t.1), silent = TRUE)

stopifnot(
          inherits(trial.1, "try-error"),
          as.logical(grep("Only one.*per side allowed", trial.1))
          )

t.2 <- t.base
t.2@trades[5,2] <- "X"

trial.2 <- try(validObject(t.2), silent = TRUE)

stopifnot(
          inherits(trial.2, "try-error"),
          as.logical(grep("Only one.*per side allowed", trial.2))
          )

## Tests for 0 shares

t.3 <- t.base
t.3@trades[1,3] <- 0

trial.3 <- try(validObject(t.3), silent = TRUE)

stopifnot(
          inherits(trial.3, "try-error"),
          as.logical(grep("Shares must be greater", trial.3))
          )

## Tests for NAs

t.4 <- t.base
t.4@trades[3,1] <- NA

trial.4 <- try(validObject(t.4), silent = TRUE)

stopifnot(
          inherits(trial.4, "try-error"),
          as.logical(grep("No NAs allowed in trades", trial.4))
          )


## Tests for invalid column names

t.5 <- t.base
names(t.5@trades)[2] <- "sides"

trial.5 <- try(validObject(t.5), silent = TRUE)

stopifnot(
          inherits(trial.5, "try-error"),
          as.logical(grep("Columns.*required", trial.5))
          )


## Tests for invalid side indicator

t.6 <- t.base
t.6@trades[1,2] <- "A"

trial.6 <- try(validObject(t.6), silent = TRUE)

stopifnot(
          inherits(trial.6, "try-error"),
          as.logical(grep("Sides must be one of", trial.6))
          )

## Tests for invalid shares class

t.7 <- t.base
t.7@trades[,3] <- as.character(t.7@trades[,3])

trial.7 <- try(validObject(t.7), silent = TRUE)

stopifnot(
          inherits(trial.7, "try-error"),
          as.logical(grep("Values.*shares.*numeric", trial.7))
          )
