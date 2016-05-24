# Tests:
library(condformat)
library(dplyr)
library(testthat)
context("show")

test_that("show_column works", {
  data(iris)
  x <- condformat(head(iris)) + show_columns(-Sepal.Length)
  expect_true("Sepal.Length" %in% colnames(x))
  out <- condformat2html(x)
  expect_that(out[1], not(matches("Sepal.Length")))
  expect_that(out[1], matches("Sepal.Width"))
  expect_that(out[1], matches("Petal.Length"))
  expect_that(out[1], matches("Petal.Width"))
  expect_that(out[1], matches("Species"))
})

test_that("show_column_ works", {
  data(iris)
  x <- condformat(head(iris)) + show_columns_(.dots = c("Sepal.Length", "Petal.Length"))
  expect_true("Sepal.Width" %in% colnames(x))
  out <- condformat2html(x)
  expect_that(out[1], not(matches("Sepal.Width")))
  expect_that(out[1], matches("Petal.Length"))
})


test_that("show_column works with custom names", {
  data(iris)
  x <- condformat(head(iris)) + show_columns(Sepal.Length, Petal.Width, Species,
                                             col_names = c("MySepLen", "MyPetWi", "MySpe"))
  expect_true("Sepal.Length" %in% colnames(x))
  expect_true("Petal.Length" %in% colnames(x))
  out <- condformat2html(x)
  expect_that(out[1], not(matches("Sepal.Length")))
  expect_that(out[1], not(matches("Sepal.Width")))
  expect_that(out[1], not(matches("Petal.Length")))
  expect_that(out[1], not(matches("Petal.Width")))
  expect_that(out[1], not(matches("Species")))
  expect_that(out[1], matches("MySepLen"))
  expect_that(out[1], matches("MyPetWi"))
  expect_that(out[1], matches("MySpe"))
})


test_that("show_row works", {
  data(iris)
  x <- condformat(head(iris, n = 10)) +
    show_rows(Sepal.Length == 5.1, Sepal.Width == 3.5,
              Petal.Length == 1.4, Petal.Width == 0.2)
  # in the data frame nothing is filtered
  expect_that(nrow(x), equals(10))
  out <- condformat2html(x)
  # the html code only shows one row (that does not have any 8 digit)
  expect_that(out[1], not(matches("8")))
})


test_that("show_row works after modifying data frame", {
  data(iris)
  x <- condformat(head(iris, n = 10))
  x$Sepal.Length <- x$Sepal.Length + 1

  x <- x + show_rows(Sepal.Length == 6.1, Sepal.Width == 3.5,
                     Petal.Length == 1.4, Petal.Width == 0.2)
  # in the data frame nothing is filtered
  expect_that(nrow(x), equals(10))
  out <- condformat2html(x)
  # the html code only shows one row (that does not have any 8 digit)
  expect_that(out[1], not(matches("8")))
})

test_that("custom show_ passes doing nothing", {
    custom_showobj <- structure(list(),
                                class = c("condformat_show_columns"))
    data(iris)
    x <- condformat(head(iris))
    y <- x + custom_showobj
    out_x <- condformat2html(x)
    out_y <- condformat2html(y)
    expect_that(out_x, is_identical_to(out_y))
})

