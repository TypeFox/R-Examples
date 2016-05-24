# Tests:
library(condformat)
library(dplyr)
library(testthat)
context("rule_fill_gradient")

test_that("rule_fill_gradient works", {
  data(iris)
  x <- condformat(iris[c(1:10, 51:60, 101:110),])
  y <- x + rule_fill_gradient(Sepal.Length)
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))

  y <- x + rule_fill_gradient_("Sepal.Length")
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))
})

test_that("rule_fill_gradient2 works", {
  data(iris)
  x <- condformat(iris[c(1:10, 51:60, 101:110),])
  y <- x + rule_fill_gradient2(Sepal.Length)
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))

  y <- x + rule_fill_gradient2_("Sepal.Length")
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))
})

