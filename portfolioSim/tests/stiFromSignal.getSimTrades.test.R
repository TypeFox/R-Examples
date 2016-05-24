################################################################################
##
## $Id: stiFromSignal.getSimTrades.test.R 413 2007-04-22 19:31:03Z enos $
##
## Tests for the getSimTrades method of class stiFromSignal
##
################################################################################

library(portfolioSim)

load("stiFromSignal.getSimTrades.test.RData")

## save(stiFS.1, stiFS.2, stiFS.3, env.1, env.2, holdings.1, holdings.2, sd.1, sd.2, sd.3, truth.1, truth.2, truth.3, truth.4, truth.5, file = "stiFromSignal.getSimTrades.test.RData", compress = TRUE)

## No holdings, no target

result.1 <- getSimTrades(stiFS.1, period = 1, holdings = new("portfolio"), sim.data = sd.1, verbose = FALSE)

stopifnot(
          all.equal(result.1@trades@trades$id, truth.1@trades@trades$id),
          all.equal(result.1@trades@trades$side, truth.1@trades@trades$side),
          all.equal(result.1@trades@trades$shares, truth.1@trades@trades$shares),
          all.equal(result.1@period, truth.1@period),
          all.equal(get("target", pos = env.1, inherits = FALSE), holdings.1)
          )

## Holdings, target, gets rebalanced

result.2 <- getSimTrades(stiFS.1, period = 1, holdings = holdings.2, sim.data = sd.1, verbose = FALSE)

stopifnot(
          all.equal(result.2@trades@trades$id, truth.2@trades@trades$id),
          all.equal(result.2@trades@trades$side, truth.2@trades@trades$side),
          all.equal(result.2@trades@trades$shares, truth.2@trades@trades$shares),
          all.equal(result.2@period, truth.2@period),
          all.equal(get("target", pos = env.1, inherits = FALSE), holdings.1)
          )

## No holdings, target, not rebalanced

result.3 <- getSimTrades(stiFS.1, period = 3, holdings = new("portfolio"), sim.data = sd.3, verbose = FALSE)

stopifnot(
          all.equal(result.3@trades@trades$id, truth.3@trades@trades$id),
          all.equal(result.3@trades@trades$side, truth.3@trades@trades$side),
          all.equal(result.3@trades@trades$shares, truth.3@trades@trades$shares),
          all.equal(result.3@period, truth.3@period)
          )

## Holdings, no target, not rebal.on (but rebalanced since no target)

result.4 <- getSimTrades(stiFS.2, period = 2, holdings = holdings.1, sim.data = sd.2, verbose = FALSE)

stopifnot(
          all.equal(result.4@trades@trades$id, truth.4@trades@trades$id),
          all.equal(result.4@trades@trades$side, truth.4@trades@trades$side),
          all.equal(result.4@trades@trades$shares, truth.4@trades@trades$shares),
          all.equal(result.4@period, truth.4@period)
          )

## Percent-volume trading style

result.5 <- getSimTrades(stiFS.3, period = 3, holdings = holdings.2, sim.data = sd.3, verbose = FALSE)

stopifnot(
          all.equal(result.5@trades@trades$id, truth.5@trades@trades$id),
          all.equal(result.5@trades@trades$side, truth.5@trades@trades$side),
          all.equal(result.5@trades@trades$shares, truth.5@trades@trades$shares),
          all.equal(result.5@period, truth.5@period)
          )

## Error: multiple periods

trial.1 <- try(
               getSimTrades(stiFS.1, period = c(1,2), holdings = holdings.2, sim.data = sd.1, verbose = FALSE),
               silent = TRUE
               )

## Error: in.var not in data

trial.2 <- try(
               getSimTrades(stiFS.2, period = 1, holdings = holdings.2, sim.data = sd.1, verbose = FALSE),
               silent = TRUE
               )

## Error: volume not in data

trial.3 <- try(
               getSimTrades(stiFS.3, period = 1, holdings = holdings.1, sim.data = sd.1, verbose = FALSE),
               silent = TRUE
               )

stopifnot(
          class(trial.1) == "try-error",
          class(trial.2) == "try-error",
          class(trial.3) == "try-error"
          )
