library(apaStyle)
context("Significance Footnote")

var1 = data.frame(
  c(8.22, 25.12, 30.27),
  c("", "+", "***")
)

test_that("Test valid dataframe", {

  expect_equal(apa.signif(data = var1)$succes, TRUE)
  expect_equal(typeof(apa.signif(data = var1)$signif), "list")

})

test_that("Test invalid dataframe", {

  test = "Invalid data is supplied."

  expect_equal(apa.signif()$succes, test)
  expect_equal(apa.signif(data = matrix(rnorm(100, mean = 0, sd = 1)))$succes, test)
  expect_equal(apa.signif(data = NULL)$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.signif(data = NULL)$signif), test)

})

test_that("Test invalid size of the dataset", {

  test = "The supplied data is too big to generate an APA formatted table."

  expect_equal(apa.signif(data = data.frame(c(1:101)))$succes, test)
  expect_equal(apa.signif(data = data.frame(matrix(runif(21), ncol = 21)))$succes, test)

  test = "NULL"

  expect_equal(typeof(apa.signif(data = data.frame(c(1:101)))$signif), test)

})
