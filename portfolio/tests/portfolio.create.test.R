################################################################################
##
## $Id: portfolio.create.test.R 1625 2010-02-18 19:44:29Z enos $
##
## Tests "create" method of "portfolio" class
##
################################################################################

library(portfolio)

load("portfolio.create.test.RData")

## save(data.0, p.0, price, truth, file = "portfolio.create.test.RData", compress = TRUE)

## Test checks that "create" method updates "shares" slot of
## "portfolio".

## erases work done by "create" method by excising data from
## the "shares" data.frame

p.0@shares <- p.0@shares[0,]

## "create" method attempts to rebuild "shares" data.frame

trial.0 <- try(
               all.equal(
                         truth$shares,
                         create(p.0)@shares$shares
                         ),
               silent = TRUE
               )

if(class(trial.0) == "try-error"){
  stop("create method failed to update shares data.frame!")
}
