library(testthat)
library(Causata)
library(ggplot2)

equals <- testthat::equals

context("BinaryPredictor")

# load diamonds data from ggplot2
data(diamonds)

result <- try(BinaryPredictor(diamonds$price, diamonds$price>5000, nbins=11, woeDelta=0, minBin=0, folds=5), silent=TRUE)

test_that("Call using arguments for BinaryCut and PredictivePowerCv works without errors.",
  expect_that(class(result), equals("BinaryPredictor")) )

test_that("Result has 11 bins in discretization.", expect_that(length(result$woe$woe.levels), equals(11)))

test_that("Result has 5 folds of cross validation.", expect_that(length(result$predictivePower$predictive.power), equals(5)))
