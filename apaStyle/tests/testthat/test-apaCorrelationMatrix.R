library(apaStyle)
context("Correlation Matrix")

var1 = data.frame(
  rnorm(100, mean = 0, sd = 1),
  rnorm(100, mean = 0, sd = 1),
  rnorm(100, mean = 0, sd = 1),
  rnorm(100, mean = 0, sd = 1)
)

test_that("Test valid dataframe", {

  expect_equal(apa.cor.matrix(data = var1)$succes, TRUE)

  expect_equal(typeof(apa.cor.matrix(data = var1)$data), "character")

  expect_equal(nrow(apa.cor.matrix(data = var1)$data), 4)
  expect_equal(ncol(apa.cor.matrix(data = var1)$data), 8)

})

test_that("Test valid position", {

  expect_equal(apa.cor.matrix(data = var1, position = "upper")$succes, TRUE)

  expect_equal(typeof(apa.cor.matrix(data = var1, position = "upper")$data), "character")

  expect_equal(nrow(apa.cor.matrix(data = var1, position = "upper")$data), 4)
  expect_equal(ncol(apa.cor.matrix(data = var1, position = "upper")$data), 8)

})

test_that("Test invalid dataframe", {

  test = "Invalid data is supplied."

  expect_equal(apa.cor.matrix(data = matrix(rnorm(100, mean = 0, sd = 1)))$succes, test)
  expect_equal(apa.cor.matrix(data = NULL)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.cor.matrix(data = NULL)$data), test)

})

test_that("Test invalid position", {

  test = "The supplied display position for the correlation matrix is not valid. Only 'upper' or 'lower' position is allowed."

  expect_equal(apa.cor.matrix(data = var1, position = "test")$succes, test)
  expect_equal(apa.cor.matrix(data = var1, position = "")$succes, test)
  expect_equal(apa.cor.matrix(data = var1, position = NULL)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.cor.matrix(data = var1, position = NULL)$data), test)

})
