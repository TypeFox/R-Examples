library(apaStyle)
context("Merge data")

var1 = rnorm(100, mean = 0, sd = 1)
var2 = rnorm(100, mean = 0, sd = 1)
var3 = header = c("M", "SD")

test_that("Test valid a values", {

  expect_equal(apa.merge(a = var2, b = var2, header = var3)$succes, TRUE)

  expect_equal(typeof(apa.merge(a = var2, b = var2, header = var3)$data), "character")
  expect_equal(typeof(apa.merge(a = var2, b = var2, header = var3)$header), "character")

  expect_equal(length(apa.merge(a = var2, b = var2, header = var3)$data), 100)
  expect_equal(length(apa.merge(a = var2, b = var2, header = var3)$header), 1)

})

test_that("Test valid y values", {

  expect_equal(apa.merge(a = var1, b = var1, header = var3)$succes, TRUE)

  expect_equal(typeof(apa.merge(a = var1, b = var1, header = var3)$data), "character")
  expect_equal(typeof(apa.merge(a = var1, b = var1, header = var3)$header), "character")

  expect_equal(length(apa.merge(a = var1, b = var1, header = var3)$data), 100)
  expect_equal(length(apa.merge(a = var1, b = var1, header = var3)$header), 1)

})

test_that("Test valid a and b values", {

  expect_equal(apa.merge(a = var1, b = var2, header = var3)$succes, TRUE)
  expect_equal(apa.merge(a = 1, b = 1, header = var3)$succes, TRUE)

  expect_equal(typeof(apa.merge(a = rnorm(100, mean = 0, sd = 1), b = rnorm(100, mean = 0, sd = 1), header = var3)$data), "character")
  expect_equal(typeof(apa.merge(a = 1, b = 1, header = var3)$data), "character")

  expect_equal(typeof(apa.merge(a = rnorm(100, mean = 0, sd = 1), b = rnorm(100, mean = 0, sd = 1), header = var3)$header), "character")
  expect_equal(typeof(apa.merge(a = 1, b = 1, header = var3)$header), "character")

  expect_equal(length(apa.merge(a = rnorm(100, mean = 0, sd = 1), b = rnorm(100, mean = 0, sd = 1), header = var3)$data), 100)
  expect_equal(length(apa.merge(a = 1, b = 1, header = var3)$data), 1)

  expect_equal(length(apa.merge(a = rnorm(100, mean = 0, sd = 1), b = rnorm(100, mean = 0, sd = 1), header = var3)$header), 1)
  expect_equal(length(apa.merge(a = 1, b = 1, header = var3)$header), 1)

})

test_that("Test invalid a values", {

  test = "Invalid data is supplied."

  expect_equal(apa.merge(a = 1, b = var2, header = var3)$succes, test)
  expect_equal(apa.merge(a = NULL, b = var2, header = var3)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.merge(a = NULL, b = var2, header = var3)$data), test)
  expect_equal(typeof(apa.merge(a = NULL, b = var2, header = var3)$header), test)

})

test_that("Test invalid b values", {

  test = "Invalid data is supplied."

  expect_equal(apa.merge(a = var1, b = 1, header = var3)$succes, test)
  expect_equal(apa.merge(a = var1, b = NULL, header = var3)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.merge(a = var1, b = NULL, header = var3)$data), test)
  expect_equal(typeof(apa.merge(a = var1, b = NULL, header = var3)$header), test)

})


test_that("Test invalid a and b values", {

  test = "Invalid data is supplied."

  expect_equal(apa.merge(a = NULL, b = NULL, header = var3)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.merge(a = NULL, b = NULL, header = var3)$data), test)
  expect_equal(typeof(apa.merge(a = NULL, b = NULL, header = var3)$header), test)

})


test_that("Test invalid headers", {

  test = "No valid headers are specified."

  expect_equal(apa.merge(a = var1, b = var2)$succes, test)
  expect_equal(apa.merge(a = var1, b = var2, header = c("var1", "var2", "var3"))$succes, test)
  expect_equal(apa.merge(a = var1, b = var2, header = NULL)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.merge(a = var1, b = var2, header = NULL)$data), test)
  expect_equal(typeof(apa.merge(a = var1, b = var2, header = NULL)$header), test)

})
