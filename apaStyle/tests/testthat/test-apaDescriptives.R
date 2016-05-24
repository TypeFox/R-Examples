library(apaStyle)
context("Generate Descriptives Table")

var1 = data.frame(
  rnorm(100, mean = 0, sd = 1),
  rnorm(100, mean = 0, sd = 1),
  rnorm(100, mean = 0, sd = 1),
  rnorm(100, mean = 0, sd = 1)
)

var2 = c("Column 1", "Column 2", "Column 3", "Column 4")

test_that("Test valid dataframe", {

  expect_equal(apa.descriptives(data = var1, variables = var2, save=FALSE)$succes, TRUE)

  expect_equal(typeof(apa.descriptives(data = var1, variables = var2, save=FALSE)$data), "list")
  expect_equal(typeof(apa.descriptives(data = var1, variables = var2, save=FALSE)$table), "list")

  expect_equal(ncol(apa.descriptives(data = var1, variables = var2, save=FALSE)$data), 11)
  expect_equal(nrow(apa.descriptives(data = var1, variables = var2, save=FALSE)$data), 4)

  expect_equal(apa.descriptives(data = var1, variables = var2, save=FALSE)$table$numcol, 11)
  expect_equal(apa.descriptives(data = var1, variables = var2, save=FALSE)$table$numrow, 4)

})

test_that("Test invalid dataframe", {

  test = "Invalid data is supplied."

  expect_equal(apa.descriptives(variables = var2)$succes, test)
  expect_equal(apa.descriptives(data = matrix(rnorm(100, mean = 0, sd = 1)), variables = var2, save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = NULL, variables = var2, save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.descriptives(data = NULL, variables = var2, save=FALSE)$data), test)
  expect_equal(typeof(apa.descriptives(data = NULL, variables = var2, save=FALSE)$table), test)

})

test_that("Test invalid variables", {

  test = "No valid variable names are specified."

  expect_equal(apa.descriptives(data = var1, save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = NULL, save=FALSE)$succes, test)

  test = "The supplied data doesn't match the specified number of variables."

  expect_equal(apa.descriptives(data = var1, variables = var2[-1], save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = append(var2, 0), save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.descriptives(data = var1, variables = NULL, save=FALSE)$data), test)
  expect_equal(typeof(apa.descriptives(data = var1, variables = NULL, save=FALSE)$table), test)

})


test_that("Test invalid items to report", {

  test = "The specified descriptives to report are not valid."

  expect_equal(apa.descriptives(data = var1, variables = var2, report = NULL, save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = var2, report = c("", ""), save=FALSE)$succes, test)

  test = "The specified descriptives to report are not valid. Only 'M', 'SD', and 'r' are allowed."

  expect_equal(apa.descriptives(data = var1, variables = var2, report = c("N", "SD", "r"), save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = var2, report = c("N", "M", "SD", "r"), save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.descriptives(data = var1, variables = var2, report = NULL, save=FALSE)$data), test)
  expect_equal(typeof(apa.descriptives(data = var1, variables = var2, report = NULL, save=FALSE)$table), test)

})

test_that("Test invalid filename", {

  test = "The supplied filename is not valid. Please specify a valid 'docx' file."

  expect_equal(apa.descriptives(data = var1, variables = var2, filename = "", save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = var2, filename = "test", save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = var2, filename = NULL, save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.descriptives(data = var1, variables = var2, filename = NULL, save=FALSE)$table), test)

})

test_that("Test invalid merge", {

  test = "The merge argument is not of logical type."

  expect_equal(apa.descriptives(data = var1, variables = var2, merge = NULL, save=FALSE)$succes, test)

  test = "Can not merge the mean and standard deviation into one column if they're not specified."

  expect_equal(apa.descriptives(data = var1, variables = var2, report = c("M"), merge = TRUE, save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = var2, report = c("SD"), merge = TRUE, save=FALSE)$succes, test)
  expect_equal(apa.descriptives(data = var1, variables = var2, report = c("M", "r"), merge = TRUE, save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.descriptives(data = var1, variables = var2, merge = NULL, save=FALSE)$table), test)

})

test_that("Test invalid size of the dataset", {

  test = "The supplied data has too many variables to generate an APA formatted table."

  expect_equal(apa.descriptives(data = data.frame(matrix(runif(21), ncol = 21)), variables = as.character(c(1:21)), save=FALSE)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.descriptives(data = data.frame(matrix(runif(21), ncol = 21)), variables = as.character(c(1:21)), save=FALSE)$table), test)

})

