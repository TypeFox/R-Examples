# Tests:
library(condformat)
library(dplyr)
library(testthat)
context("rule_fill_discrete")

test_that("rule_fill_discrete works", {
  data(iris)
  x <- condformat(iris[c(1:10, 51:60, 101:110),])
  y <- x + rule_fill_discrete(Species,
                              expression  = Sepal.Length > max(Sepal.Length),
                              colours = c("TRUE" = "red", "FALSE" = "blue"))
  out <- condformat2html(y)
  expect_that(out[1], not(matches("red")))

  y <- x + rule_fill_discrete(Species,
                              expression  = Sepal.Length >= min(Sepal.Length),
                              colours = c("TRUE" = "red", "FALSE" = "blue"))
  out <- condformat2html(y)
  expect_that(out[1], not(matches("blue")))

  y <- x + rule_fill_discrete(Species)
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))
})


test_that("rule_fill_discrete lock cells", {
  data(iris)
  x <- condformat(head(iris))
  y <- x +
    rule_fill_discrete(Species,
                       expression  = 1,
                       colours = c("1" = "red")) +
    rule_fill_discrete(Species,
                       expression  = 1,
                       colours = c("1" = "blue"))
  out <- condformat2html(y)
  expect_that(out[1], not(matches("red")))

  y <- x +
    rule_fill_discrete(Species,
                       expression  = 1,
                       colours = c("1" = "red"),
                       lockcells = TRUE) +
    rule_fill_discrete(Species,
                       expression  = 1,
                       colours = c("1" = "blue"))
  out <- condformat2html(y)
  expect_that(out[1], not(matches("blue")))

})


test_that("rule_fill_discrete(_) syntax with multiple variables and no expression gives warning", {
  expect_that(rule_fill_discrete(Species, Sepal.Length),
              gives_warning("multiple variables"))
  expect_that(rule_fill_discrete_(columns = c("Species", "Sepal.Length")),
              gives_warning("multiple variables"))
})

test_that("rule_fill_discrete_ works", {
  data(iris)
  x <- condformat(head(iris))
  y <- x + rule_fill_discrete_("Species")
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))
  y <- x + rule_fill_discrete_("Species", expression = "Sepal.Length > 4.6",
                               colours = c("TRUE" = "red"))
  out <- condformat2html(y)
  expect_that(out[1], matches("^<table.*</table>$"))
})

test_that("custom rule_ passes doing nothing", {
  custom_ruleobj <- structure(list(),
                              class = c("condformat_rule"))
  data(iris)
  x <- condformat(head(iris))
  y <- x + custom_ruleobj
  out_x <- condformat2html(x)
  out_y <- condformat2html(y)
  expect_that(out_x, is_identical_to(out_y))
})
