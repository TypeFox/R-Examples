################################################################################
##
## $Id: portfolio.expose.test.R 406 2007-04-19 16:30:22Z enos $
##
## Tests the expose method of class portfolio.
##
################################################################################

library(portfolio)

load("portfolio.expose.test.RData")

## save(test.portfolio, test.trades, empty.trades, bad.1, bad.2, bad.3, bad.4, truth, file = "portfolio.expose.test.RData", compress = TRUE)

result.1 <- expose(test.portfolio, test.trades)
result.1 <- calcWeights(result.1)

result.2 <- expose(test.portfolio, empty.trades)
result.2 <- calcWeights(result.2)

stopifnot(
          all.equal(result.1, truth),
          all.equal(result.2, test.portfolio)
          )


## Try to sell when you are shorting

trial.1 <- try(
               expose(test.portfolio, bad.1), silent = TRUE
               )

stopifnot(
          inherits(trial.1, "try-error"),
          as.logical(grep("Illegal trades found", trial.1[1]))
          )

## Try to short when you are already longing

trial.2 <- try(
               expose(test.portfolio, bad.2), silent = TRUE
               )

stopifnot(
          inherits(trial.2, "try-error"),
          as.logical(grep("Illegal trades found", trial.2[1]))
          )

## Try to cover more than you are shorting

trial.3 <- try(
               expose(test.portfolio, bad.3), silent = TRUE
               )

stopifnot(
          inherits(trial.3, "try-error"),
          as.logical(grep("Illegal trades found", trial.3[1]))
          )

## Try a side change, from long to short, without selling all your holdings

trial.4 <- try(
               expose(test.portfolio, bad.4), silent = TRUE
               )

stopifnot(
          inherits(trial.4, "try-error"),
          as.logical(grep("Illegal trades found", trial.4[1]))
          )
