################################################################################
##
## $Id: classes.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests explicitly declared validity functions in AllClasses.R
##
################################################################################

library(portfolio)

## constructs lists that break the validity functions of objects of
## class "objectHistory" and objects of class "portfolio"

list.0 <- list(data.frame(), new("portfolio"))
list.1 <- list(new("contribution"), new("contribution"))
list.2 <- list(new("performance"), new("performance"))
list.3 <- list(new("exposure"), new("exposure"))

## constucts a data.frame to test the validity function of "portfolio"

data.0 <- data.frame(id = 1:20, in.var = 1:20, ret.var = 1:20, price.var = 1:20)
p <- new("portfolio", data = data.0, id.var = "id",
              in.var = "in.var", price.var = "price.var",
              ret.var = "ret.var", sides = c("long", "short"), equity = 10000)

## directly modify "shares" data.frame so validity function returns
## false

p@shares <- p@shares[1,]

## tests validity functions of forementioned classes

trial.0 <- try(
               new("objectHistory", data = list.0), silent = TRUE
               )

trial.1 <- try(
               new("performanceHistory", data = list.1), silent = TRUE
               )

trial.2 <- try(
               new("exposureHistory", data = list.2), silent = TRUE
               )

trial.3 <- try(
               new("contributionHistory", data = list.3), silent = TRUE
               )

trial.4 <- try(
               validObject(p), silent = TRUE
               )

if(class(trial.0) == "try-error"){
  stopifnot(
            as.logical(grep("Error.*validObject.*objectHistory",trial.0[1]))
            )
}

if(class(trial.1) == "try-error"){
  stopifnot(
            as.logical(grep("Error.*validObject.*performanceHistory",
                            trial.1[1]))
            )
}

if(class(trial.2) == "try-error"){
  stopifnot(
            as.logical(grep("Error.*validObject.*exposureHistory",trial.2[1]))
            )
}

if(class(trial.3) == "try-error"){
  stopifnot(
            as.logical(grep("Error.*validObject.*contributionHistory",
                            trial.3[1]))
            )
}

if(class(trial.4) == "try-error"){
  stopifnot(
            as.logical(grep("Error.*validObject.*portfolio",
                trial.4[1]))
            )
}

## Tests the 'validity' method of 'matchedPortfolio'

test <- try(
            new("matchedPortfolio", formula = y ~ x + z, original = p),
            silent = TRUE
            )

stopifnot(
          as.logical(grep("Error.*validObject.*does not contain columns",
                          test[1]))
          )
