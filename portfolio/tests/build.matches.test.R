################################################################################
##
## $Id: build.matches.test.R 346 2006-10-01 05:08:55Z enos $
##
## tests the .build.matches private function
##
################################################################################

library(portfolio)

universe  <- c("IBM", "AAPL", "MSFT", "GOOG", "DELL", "HPQ", "PFE",
              "MRK", "DOW", "HAL")

id.matrix <- matrix(rep("HPQ", 25),
                    nrow = 5,
                    ncol = 5,
                    dimnames = list(universe[1:5], 1:5)
                    )

weights <- rep(0.1, 5)

matches <- portfolio:::.build.matches(universe, id.matrix, weights)

truth <- matrix(0,
                nrow = 10,
                ncol = 5,
                dimnames = list(universe, 1:5)
                )

truth["HPQ",] <- 0.5
  
stopifnot(
          all.equal(matches, truth)
          )
