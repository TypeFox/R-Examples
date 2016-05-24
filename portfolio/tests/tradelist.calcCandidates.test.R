################################################################################
##
## $Id: tradelist.calcCandidates.test.R 346 2006-10-01 05:08:55Z enos $
##
## Tests calcCandidates method of "tradelist" class
##
################################################################################

library(portfolio)

load("tradelist.calcCandidates.test.RData")

## save(data, p.orig, p.target, truth.candidates, truth.restricted, truth.2, truth.3, file = "tradelist.calcCandidates.test.RData", compress = TRUE)

## Creates base objects for use in the remainder of the test

tl.base <- new("tradelist", turnover = 64000, data = data)

## Simulates actual use of the tradelist class by using 2 non-empty
## portfolio.  Tests the overall functionality of the method.

tl.1 <- portfolio:::calcCandidates(tl.base, p.orig, p.target)

stopifnot(all.equal(tl.1@candidates, truth.candidates))
stopifnot(all.equal(tl.1@restricted, truth.restricted))

## Tests appropriate behavior for when one the "original" portfolio is
## empty and the target portfolio contains positions

p.empty <- new("portfolio")

tl.2 <- portfolio:::calcCandidates(tl.base, p.empty, p.target)

## Tests the corner case where both portfolios are empty

trial.1 <- try(
               new("tradelist", p.empty, p.empty,
                   turnover = 1000, data = data),
               silent = TRUE
               )

stopifnot(class(trial.1) == "try-error" && 
          as.logical(grep("Both portfolios.*empty", trial.1[1]))
)

## Tests the overall behavior of a tradelist when the "unrestricted"
## flag has been set and type has been set to "all"

tl.base.unrestricted <- new("tradelist", turnover = 64000, data =
                            data, type = "all", unrestricted = TRUE)

tl.3 <- portfolio:::calcCandidates(tl.base.unrestricted, p.orig, p.target)
tl.3@candidates <- tl.3@candidates[order(row.names(tl.3@candidates)),]

stopifnot(all.equal(tl.2@candidates$shares, truth.2$shares))
stopifnot(all.equal(tl.2@candidates$side, truth.2$side))
stopifnot(all.equal(tl.2@candidates$mv, truth.2$mv))

## Sets the id.var slot to "id.var". Proper behavior includes creating
## another column in data named "id" with the same values as the
## "id.var" column

tl.4 <- tl.base

tl.4@id.var <- "id.var"
names(tl.4@data)[match("id", names(tl.4@data))] <- "id.var"

tl.4 <- portfolio:::calcCandidates(tl.4, p.orig, p.target)

stopifnot(all.equal(tl.4@data$id.var, tl.4@data$id))

## If a stock is in the "shares" slot of one or both of the
## portfolios, but there is not a row for it in the "data" slot of the
## "tradelist" object, that stock should be put on the restricted
## list.

## Create a "tradelist", removing certain rows from "data"

tl.5 <- new("tradelist", turnover = 10000, data = data[-c(1:3,6,16),], unrestricted = FALSE)

tl.5 <- portfolio:::calcCandidates(tl.5, p.orig, p.target)

## if restrictions work properly, the stocks which had their rows
## removed from "data" will be in the "restricted" data frame and no
## longer appear in the "candidates" data frame

stopifnot(
          all(c(2,3) %in% tl.5@restricted$id),
          !any(c(2,3) %in% tl.5@candidates$id)
          )

## If there is missing data (NAs) in data for a stock, that stock
## should be added to the restricted list

tl.6 <- tl.base

## sets the prices of certain stocks to NA

tl.6@data$price.usd[2:5] <- NA

tl.6 <- portfolio:::calcCandidates(tl.6, p.orig, p.target)

## Stocks "2:5" are all candidates, but should be placed on the
## restriced list because there is no price data for them

stopifnot(
          all(c(2:5) %in% tl.6@restricted$id),
          !any(c(2:5) %in% tl.6@candidates$id)
          )

## If there is missing data (volume in this case) in data for a stock,
## that stock should be added to the restricted list

tl.7 <- tl.base
tl.7@data$volume[c(2:5)] <- NA

tl.7 <- portfolio:::calcCandidates(tl.7, p.orig, p.target)

## Stocks "2:5" are all candidates, but should be placed on the
## restriced list because there is no volume data for them

stopifnot(
          all(c(2:5) %in% tl.7@restricted$id),
          !any(c(2:5) %in% tl.7@candidates$id)
          )

## tests restrictions on trades with shares below trade.usd.min slot

tl.8 <- tl.base
tl.8@trade.usd.min <- 4500

tl.8 <- portfolio:::calcCandidates(tl.8, p.orig, p.target)

stopifnot(
          all(c(2,3) %in% tl.8@restricted$id),
          !any(c(2,3) %in% tl.8@candidates$id)
          )



stopifnot(all.equal(row.names(tl.3@candidates), truth.3))

## Construction of data frame to be used throughout tradelist tests.
## We include these steps as reference for how we created the initial
## set of tests.  The data frame stored in the data slot in objects in
## subsequent tests may be different than this initial version. We
## build portfolios so that on the long and short sides 1 position
## does not change, 1 position decreases, 1 position increases, 1
## position closes, 1 position opens, and 1 position changes sides.


## data              <- data.frame(id = 1:20)
## data$volume       <- seq(2100, 4000, by = 100)
## data$price.usd    <- seq(10, 200, by = 10)
## data$round.lot    <- 1
## data$sort.input.1 <- rnorm(1:nrow(data))
## data$sort.input.2 <- rnorm(1:nrow(data), mean = 1, sd = 4)

## data$p.orig.shares   <- NA
## data$p.target.shares <- NA

## data[data$id == 1,c("p.orig.shares","p.target.shares")] <- c(100, 100)
## data[data$id == 2,c("p.orig.shares","p.target.shares")] <- c(100, 50)
## data[data$id == 3,c("p.orig.shares","p.target.shares")] <- c(100, 150)
## data[data$id == 4,c("p.orig.shares","p.target.shares")] <- c(100, NA)
## data[data$id == 5,c("p.orig.shares","p.target.shares")] <- c(NA, 100)
## data[data$id == 6,c("p.orig.shares","p.target.shares")] <- c(-50, 100)

## data[data$id == 11,c("p.orig.shares","p.target.shares")] <- c(-100, -100)
## data[data$id == 12,c("p.orig.shares","p.target.shares")] <- c(-100, -50)
## data[data$id == 13,c("p.orig.shares","p.target.shares")] <- c(-100, -150)
## data[data$id == 14,c("p.orig.shares","p.target.shares")] <- c(-100, NA)
## data[data$id == 15,c("p.orig.shares","p.target.shares")] <- c(NA, -100)
## data[data$id == 16,c("p.orig.shares","p.target.shares")] <- c(50, -100)

## data$id <- as.character(data$id)

## p.orig <- new("portfolio", data = data, price.var = "price.usd")
## shares <- data[!is.na(data$p.orig.shares), c("id","p.orig.shares")]
## names(shares) <- c("id","shares")
## p.orig@shares <- shares
## p.orig <- calcWeights(p.orig)

## p.target <- new("portfolio", data = data, price.var = "price.usd")
## shares <- data[!is.na(data$p.target.shares),c("id","p.target.shares")]
## names(shares) <- c("id","shares")
## p.target@shares <- shares
## p.target <- calcWeights(p.target)

## Miscellaneous Useful Data

## mv.orig <- sum(tl@data$price.usd * abs(tl@data$p.orig.shares), na.rm = TRUE)
## mv.target <- sum(tl@data$price.usd * abs(tl@data$p.target.shares), na.rm = TRUE)

