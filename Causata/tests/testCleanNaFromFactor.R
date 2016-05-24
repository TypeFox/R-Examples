# Replaces missing values in a factor
# 
# Author: Justin

library(testthat)
library(Causata)

context("CleanNaFromFactor")

# create a factor with missing values
f <- as.factor(c("a","b","c",NA))
test_that("Missing values are present in test data", expect_true( any(is.na(f)) ))
test_that("Missing values are replaced", expect_false( any(is.na(CleanNaFromFactor(f))) ))