library(apaStyle)
context("Generate Table")

var1 = data.frame(
  c("Column 1", "Column 2", "Column 3"),
  c(3.45, 5.21, 2.64),
  c(1.23, 1.06, 1.12),
  c(8.22, 25.12, 30.27),
  c("+", "**", "***")
)

var2 = c("Variable", "M", "SD", "t-value", "*")

test_that("Test valid dataframe", {

  expect_equal(apa.table(data = var1, level1.header = var2, save=FALSE)$succes, TRUE)

  expect_equal(typeof(apa.table(data = var1, level1.header = var2, save=FALSE)$table), "list")

  expect_equal(apa.table(data = var1, level1.header = var2, save=FALSE)$table$numcol, 5)
  expect_equal(apa.table(data = var1, level1.header = var2, save=FALSE)$table$numrow, 3)

})

test_that("Test invalid dataframe", {

  test = "Invalid data is supplied."

  expect_equal(apa.table(level1.header = var2)$succes, test)
  expect_equal(apa.table(data = matrix(rnorm(100, mean = 0, sd = 1)), level1.header = var2, save=FALSE)$succes, test)
  expect_equal(apa.table(data = NULL, level1.header = var2, save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.table(data = NULL, level1.header = var2, save=FALSE)$table), test)

})

test_that("Test invalid header", {

  test = "No valid headers are specified."

  expect_equal(apa.table(data = var1, save=FALSE)$succes, test)
  expect_equal(apa.table(data = var1, level1.header = NULL, save=FALSE)$succes, test)

  test = "The level 1 colspan doesn't match the number of level 2 headers."

  expect_equal(apa.table(data = var1, level1.header = c("Descriptives", "Inferential"), level2.header = var2, save=FALSE)$succes, test)
  expect_equal(apa.table(data = var1, level1.header = c("Descriptives", "Inferential"), level1.colspan = c(1, 1), level2.header = var2, save=FALSE)$succes, test)
  test = "NULL"

  expect_equal(typeof(apa.table(data = var1, level1.header = NULL, save=FALSE)$table), test)

})

test_that("Test invalid filename", {

  test = "The supplied filename is not valid. Please specify a valid 'docx' file."

  expect_equal(apa.table(data = var1, level1.header = var2, filename = "", save=TRUE)$succes, test)
  expect_equal(apa.table(data = var1, level1.header = var2, filename = "test", save=TRUE)$succes, test)
  expect_equal(apa.table(data = var1, level1.header = var2, filename = NULL, save=TRUE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.table(data = var1, level1.header = var2, filename = NULL, save=TRUE)$table), test)

})

test_that("Test invalid size of the dataset", {

  test = "The supplied data is too big to generate an APA formatted table."

  expect_equal(apa.table(data = data.frame(c(1:101)), level1.header = as.character(c(1)), save=FALSE)$succes, test)
  expect_equal(apa.table(data = data.frame(matrix(runif(21), ncol = 21)), level1.header = as.character(c(1:21)), save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.table(data = data.frame(c(1:101)), level1.header = as.character(c(1:20)), save=FALSE)$table), test)

})

test_that("Test if the length of the dataframe does mot match with the length of the header", {

  test = "The supplied data doesn't match the specified table header."

  expect_equal(apa.table(data = var1, level1.header = var2[-1], save=FALSE)$succes, test)
  expect_equal(apa.table(data = var1, level1.header = append(var2, 0), save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.table(data = var1, level1.header = var2[-1], save=FALSE)$table), test)

})


