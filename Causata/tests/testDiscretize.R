library(testthat)
library(Causata)

#source("utils.R")
equals <- testthat::equals

context("Discretize")

x.nomissing <- c(1.1,1.2,1.3, 2.1,2.2,2.3, 3.1,3.2,3.3, 4.1,4.2,4.3)
x.missing   <- c(1.1,1.2,1.3, 2.1,2.2,2.3, 3.1,3.2,3.3,  NA, NA, NA)
breaks    <- c(1,2,3,4,5)
breaks.na <- c(1,2,3,4)
discrete.values    <- c(10,20,30,40)
discrete.values.na <- c(10,20,30, 0)
causataData <- CausataData(data.frame('x__AP'=x.nomissing, xna__AP=x.missing), rep(0,length(x.nomissing)))

test_that("Not using ReplaceOutliers throws error",
  expect_that(Discretize(causataData, 'x__AP', breaks, discrete.values), throws_error() )
)

causataData2 <- ReplaceOutliers(causataData, 'x__AP', lowerLimit=1)
test_that("Setting lower limit with ReplaceOutliers but not upper limit throws error",
  expect_that(Discretize(causataData2, 'x__AP', breaks, discrete.values), throws_error() )
)

causataData3 <- CleanNaFromContinuous(causataData, 'x__AP')
test_that("Cleaning missing values before discretization throws error",
  expect_that(Discretize(causataData3, 'x__AP', breaks, discrete.values), throws_error() )
)

causataData4 <- ReplaceOutliers(causataData, 'x__AP', lowerLimit=min(causataData$df$x__AP), upperLimit=max(causataData$df$x__AP))
causataData4 <- Discretize(causataData4, 'x__AP', breaks, discrete.values)
test_that("Values are mapped to discrete values",
  expect_equal(causataData4$df$x__AP, c(10,10,10, 20,20,20, 30,30,30, 40,40,40))
)

causataData5 <- ReplaceOutliers(causataData, 'xna__AP', lowerLimit=min(causataData$df$xna__AP, na.rm=TRUE), upperLimit=max(causataData$df$xna__AP, na.rm=TRUE))
causataData5 <- Discretize(causataData5, 'xna__AP', breaks.na, discrete.values.na)
test_that("Values are mapped to discrete values, missing mapped to zero",
  expect_equal(causataData5$df$xna__AP, c(10,10,10, 20,20,20, 30,30,30, 0,0,0))
)