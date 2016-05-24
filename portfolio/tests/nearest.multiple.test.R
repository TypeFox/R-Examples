################################################################################
##
## $Id: nearest.multiple.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests for the "nearest.multiple" function
##
################################################################################

library(portfolio)

load("nearest.multiple.test.RData")

## save(truth, file = "nearest.multiple.test.RData", compress = TRUE)

## Invalid input: complex numbers

trial.0 <- try(
               portfolio:::.nearest.multiple(5+0i, 10), silent = TRUE
               )

if(class(trial.0) == "try-error"){
  stopifnot(
            as.logical(grep("Error.*invalid.*complex",trial.0[1]))
            )
}

stopifnot(
          all.equal(truth$result,
                    portfolio:::.nearest.multiple(truth$x, truth$y))
          )
