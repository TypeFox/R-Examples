# Test variable transformation replay
# Author: David Barker
###############################################################################
library(testthat)
library(Causata)

source("utils.R")
equals <- testthat::equals

context("TransformationReplay")

test_that("Applying MergeLevels can be replayed", {
  max.levels <- 3
  
  df <- data.frame(var__AP=c("a", "a", "a",  "b", "b",  "c", "d", "e", NA))
  cdata <- CausataData(df, rep(0,nrow(df)))
  cdata <- MergeLevels(cdata, max.levels=max.levels)
  transformer <- GetTransforms(cdata)
  
  # First, let's be sure that the merging has done what we expect.
  #expect_that(cdata$df$var__AP, is_identical_to(factor(c("a", "a", "a",  "b", "b",  "Other", "Other", "Other", NA))))
  f1 <- cdata$df$var__AP
  f2 <- factor(c("a", "a", "a",  "b", "b",  "Other", "Other", "Other", NA))
  expect_that(f1, is_equivalent_to(f2))
  
  # Now, create a new data.frame and apply the same transformation to it.  Any values in the factor we have
  # that is not "a", or "b" should be replaced with "Other"
  #
  new.df <- data.frame(var__AP=c("a", "b", "c", "c", "c", "d", "a", NA))
  transformed <- transformer(new.df)
  expect_that(transformed$var__AP, is_equivalent_to(factor(c("a", "b", "Other", "Other", "Other", "Other", "a", NA))))
})

test_that("Applying CleanNaFromFactor can be replayed", {
  df <- data.frame(var__AP=c("a", "a", "a", NA,  "b", "b",  NA, "c", "d", "e", NA))
  cdata <- CausataData(df, rep(0,nrow(df)))
  cdata <- CleanNaFromFactor(cdata, replacement="RS")
  transformer <- GetTransforms(cdata)
  
  expect_that(cdata$df$var__AP, is_equivalent_to(factor(c("a", "a", "a", "RS", "b", "b", "RS", "c", "d", "e", "RS"))))
  
  new.df <- data.frame(var__AP=c("a", "a", NA, "b", NA, NA))
  expected <- factor(c("a", "a", "RS", "b", "RS", "RS"))
  actual <- transformer(new.df)$var__AP
  
  expect_that(actual, is_equivalent_to(expected))
})

test_that("Applying CleanNaFromContinuous can be replayed", {
  input <- c(1, 1.2, NA,  2, 3,  NA, 4, 4.5, NA)
  df <- data.frame(var__AP=input)
  cdata <- CausataData(df, rep(0,nrow(df)))
  cdata <- CleanNaFromContinuous(cdata, method="mean")
  transformer <- GetTransforms(cdata)
  
  input.mean <- mean(input, na.rm=TRUE)
  expect_that(cdata$df$var__AP, equals(c(1, 1.2, input.mean, 2, 3, input.mean, 4, 4.5, input.mean)))
  
  new.df <- data.frame(var__AP=c(1, NA, 10, NA, 20, NA))
  expected <- c(1, input.mean, 10, input.mean, 20, input.mean)
  actual <- transformer(new.df)$var__AP
  
  expect_that(actual, equals(expected))
})

test_that("Applying ReplaceOutliers can be replayed", {
  input <- c(-3, -2, -1, 0, 1, 2, 3)
  df <- data.frame(var__AP=input)
  cdata <- CausataData(df, rep(0,nrow(df)))
  cdata <- ReplaceOutliers(cdata, "var__AP", lowerLimit=-2, upperLimit=2)
  transformer <- GetTransforms(cdata)

  expect_that(cdata$df$var__AP, equals(c(-2, -2, -1, 0, 1, 2, 2)))
  
  new.df <- data.frame(var__AP=c(-4, -3, -2, 1, 2, 3, 4))
  expected <- c(-2, -2, -2, 1, 2, 2, 2)
  actual <- transformer(new.df)$var__AP
  
  expect_that(actual, equals(expected))
})

test_that("Two transformations to one column can be replayed", {
  # Two transforms will be applied
  # first replace missing values with the median value
  # then replace outliers
  #
  input.column <- c(NA, 1:10, NA, 11:21, NA)
  df <- data.frame(var__AP=input.column)
  cdata <- CausataData(df, rep(0,nrow(df)))
  
  # Transform 1
  cdata <- CleanNaFromContinuous(cdata, method="median")
  expect_that(sum(is.na(cdata$df$var__AP)), equals(0))
  
  # Transform 2
  cdata <- ReplaceOutliers(cdata, 'var__AP', upperLimit=19)
  expect_that(max(cdata$df$var__AP), equals(19))
  
  transformer <- GetTransforms(cdata)
  
  # Now create the same data frame again (with the same input) and expect that the transformation funciton
  df.2 <- data.frame(var__AP=input.column, dependent.variable=rep(0,length(input.column)))
  df.2.transformed <- transformer(df.2)
  expect_that(df.2.transformed$var__AP, equals(cdata$df$var__AP))
})

test_that("Two transformations to two columns can be replayed", {
  input.column.1 <- c(NA, 1:10, 10, 11:21)
  input.column.2 <- c(NA, 1:10, 11, 11:21)
  df <- data.frame(var1__AP=input.column.1, var2__AP=input.column.2)
  cdata <- CausataData(df, rep(0,nrow(df)))

  # don't specify a variable, so this is applied to all numeric columns
  cdata <- CleanNaFromContinuous(cdata)
  expect_that(sum(is.na(cdata$df$var1__AP)), equals(0))
  expect_that(sum(is.na(cdata$df$var2__AP)), equals(0))
  
  transformer <- GetTransforms(cdata)
  
  # Now create the same data frame again (with the same input) and expect that the transformation funciton
  df.2 <- data.frame(var1__AP=input.column.1, var2__AP=input.column.2, dependent.variable=rep(0,length(input.column.1)))
  df.2.transformed <- transformer(df.2)
  expect_that(df.2.transformed$var1__AP, equals(cdata$df$var1__AP))
  expect_that(df.2.transformed$var2__AP, equals(cdata$df$var2__AP))
})

test_that("Applying Discretize can be replayed", {
  x <- c(1.1,1.2,1.3, 2.1,2.2,2.3, 3.1,3.2,3.3,  NA, NA, NA)
  breaks <- c(1,2,3,4)
  discrete.values <- c(10,20,30,0)
  df <- data.frame(x__AP=x, dependent.variable=rep(0,length(x)))
  causataData <- CausataData(df)
  causataData <- ReplaceOutliers(causataData, 'x__AP', lowerLimit=min(x, na.rm=TRUE), upperLimit=max(x, na.rm=TRUE))
  causataData <- Discretize(causataData, 'x__AP', breaks, discrete.values)
  
  transformer <- GetTransforms(causataData)
  
  # apply the transformer function to two data frames
  # - the original -- results should match
  # - a variation on the original where values above the outlier limit are provided
  #   this tests if the outlier is mapped to the last non-missing discrete value (30)
  df.new <- data.frame(
    x__AP    = c(1.1,1.2,1.3, 2.1,2.2,2.3, 3.1,5,100,  NA, NA, NA),
    expected = c(10 ,10 , 10, 20 ,20 ,20 , 30 ,30,30,   0,  0,  0),
    dependent.variable=rep(0, length(x)))
  df.transformed     <- transformer(df)
  df.transformed.new <- transformer(df.new)
  
  # confirm that values are mapped correctly
  expect_equal(causataData$df$x__AP, df.transformed$x__AP, label="Transformed data matches original")
  expect_equal(df.transformed.new$x__AP, df.new$expected, label="New transformed data handles outliers correctly")
})