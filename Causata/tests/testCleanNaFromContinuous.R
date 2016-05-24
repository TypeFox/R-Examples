# Replaces missing values in vectors
# 
# Author: Justin

library(testthat)
library(Causata)

context("CleanNaFromContinuous")

# create a data frame with missing values
var1 <- as.integer(c(NA, 1,1  ,2,3, 8)) # integer, median is 2, mean is 3
var2 <-            c(NA, 1,1.5,2,3,11)  # numeric, median is 2, mean is 3.7
var1.list <- CleanNaFromContinuous(var1, return.replacement=TRUE)
var2.list <- CleanNaFromContinuous(var2, method="mean", return.replacement=TRUE)

test_that("Default method replaces NA with median values", expect_equal( 2, var1.list$replacement.value ))
test_that("Mean method replaces NA with mean values", expect_equal( 3.7, var2.list$replacement.value ))
test_that("Sending bad method produces errors", expect_that( CleanNaFromContinuousVariables(df, method="qwerty"), throws_error() ))
