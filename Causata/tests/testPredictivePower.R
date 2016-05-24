
library(testthat)
library(Causata)
library(stringr)

context("PredictivePower")

equals <- testthat::equals

test_that("Power is 1/3 where levels differ by 1/3, missing values in iv are ignored.",
          expect_that(PredictivePower(
            factor(c(str_split("a a a b b b", " ")[[1]], NA,NA)),
            c(                  1,1,0,0,0,1,              1, 1 ) ), equals(1/3))
)

test_that("Power is 1.0 for perfect predictor.",
          expect_that(PredictivePower(
            factor(c(str_split("a a a a a b b b b b", " "))[[1]]),
            factor(c(str_split("1 1 1 1 1 0 0 0 0 0", " "))[[1]]) ), equals(1))
)

test_that("Power is 0 for random predictor.",
          expect_that(PredictivePower(
            factor(c(str_split("a a a a b b b b", " "))[[1]]),
            factor(c(str_split("1 1 0 0 1 1 0 0", " "))[[1]]) ), equals(0))
)

test_that("Mismatched lengths throw error.",
          expect_that(PredictivePower(
            factor(c(str_split("a a a a b b b b", " "))[[1]]),
            factor(c(str_split("1 1 0 0 1 1 0",   " "))[[1]]) ), throws_error())
)

test_that("Missing value in DV throws error.",
          expect_that(PredictivePower(
            factor(c(str_split("a a b b", " "))[[1]]),
            c(NA,1,0,1) ), throws_error())
)

test_that("Number of classes not equal to two in DV throws error.", {
          expect_that(PredictivePower(
            factor(c(str_split("a a b b", " "))[[1]]),
            c(1,2,3,4) ), throws_error()) # four classes
          expect_that(PredictivePower(
            factor(c(str_split("a a b b", " "))[[1]]),
            c(1,1,1,1) ), throws_error()) } # one class
)
