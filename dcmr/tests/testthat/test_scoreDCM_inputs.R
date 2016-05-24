context('ScoreDCM inputs')

test_that("Test ScoreDCM gives error for missing arguments", {
  expect_error(ScoreDCM(), "argument \"observations\" is missing, with no default")
  expect_error(ScoreDCM(observations = 1), "argument \"qmatrix\" is missing, with no default")
  expect_error(ScoreDCM(observations = 1, qmatrix = 1), "argument \"parameter.means\" is missing, with no default")
})

test_that('Test ScoreDCM gives error for NULL observations', {
  expect_error(ScoreDCM(observations = NULL, qmatrix = 1, parameter.means = 2), "argument \"observations\" must be data frame or a matrix")
  expect_error(ScoreDCM(observations = "", qmatrix = 1, parameter.means = 2), "argument \"observations\" must be data frame or a matrix")
})

test_that('Test ScoreDCM gives error for empty observations', {
  expect_error(ScoreDCM(observations = data.frame(), qmatrix = 1, parameter.means = 2), "argument \"observations\" can not be empty")
  expect_error(ScoreDCM(observations = matrix(), qmatrix = 1, parameter.means = 2))
})

test_that('Test ScoreDCM throws error for non-numeric and non binary value of observations', {
  expect_error(ScoreDCM(observations = data.frame("hello"), qmatrix = 1, parameter.means = 2),
               "argument \"observations\" must be a matrix of all numeric values")
  expect_error(ScoreDCM(observations = matrix("hello"), qmatrix = 1, parameter.means = 2),
               "argument \"observations\" must be a matrix of all numeric values")
  expect_error(ScoreDCM(observations = data.frame(2), qmatrix = 1, parameter.means = 2),
               "argument \"observations\" must be a matrix of 0s and 1s")
  expect_error(ScoreDCM(observations = matrix(2), qmatrix = 1, parameter.means = 2),
               "argument \"observations\" must be a matrix of 0s and 1s")
})

test_that('Test ScoreDCM throws error for empty, non-numeric and non binary value of qmatrix',{
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = NULL, parameter.means = 2))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(), parameter.means = 2))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix("hello"), parameter.means = 2))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = data.frame(), parameter.means = 2))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(4), parameter.means = 2))
})

test_that('Test ScoreDCM throws error for incorrect size of qmatrix', {
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(1), parameter.means = 2))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(c(1, 1), byrow=FALSE), parameter.means = 2))
})

test_that('Test ScoreDCM throws error for empty or non numeric value of parameter.means', {
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(1), parameter.means = "hello"))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(1), parameter.means = ""))
  expect_error(ScoreDCM(observations = matrix(c(0, 1)), qmatrix = matrix(1), parameter.means = NULL))
})





