
library(testthat)
library(Causata)
library(stringr)

equals <- testthat::equals

context("Woe")

# create a factor with three levels
# - odds of 1 for a:  2:1 = 2.0
# - odds of 1 for b:  1:2 = 0.5
# - odds of 1 for NA: 1:1 = 1.0
f  <- factor(c(str_split("a a a b b b", " ")[[1]], NA,NA))
dv <- c(                  1,1,0,0,0,1,              1, 0 )
fw <- Woe(f,dv)


test_that("Odds and log odds match expected values, and missing values are handled.", {
          expect_that(fw$odds, equals(c(2.0,0.5,1.0)))
          expect_equivalent(fw$woe.levels, log(c(2.0,0.5,1.0))) }
)
