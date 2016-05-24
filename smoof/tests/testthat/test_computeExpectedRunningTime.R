context("Expected Running Time")

test_that("computeExpectedRunningTime computes ERT the right way.", {
  # check that error is thrown on wrong parameter settings
  expect_error(computeExpectedRunningTime(
    fun.evals = c(100L, 1000L, 10000L),
    fun.success.runs = c(TRUE, TRUE, FALSE),
    fun.reached.target.values = c(0.1, 0.01, 0.5)
  ))

  # missing fun.target.value
  expect_error(computeExpectedRunningTime(
    fun.evals = c(100L, 200L, 600L),
    fun.reached.target.values = c(0.01, 0.001, 0.075)
  ))

  # check that values are correct if we pass evals and success vector
  ERT = computeExpectedRunningTime(
    fun.evals = c(100L, 200L, 600L),
    fun.success.runs = c(TRUE, TRUE, TRUE)
  )
  expect_true(is.numeric(ERT))
  expect_equal(ERT, 300L)

  # check that vector of success is computed internally
  ERT = computeExpectedRunningTime(
    fun.evals = c(100L, 200L, 600L),
    fun.reached.target.values = c(0.01, 0.001, 0.075),
    fun.target.value = 0.1
  )
  expect_equal(ERT, 300L)

  # check penalty value
  ERT = computeExpectedRunningTime(
    fun.evals = rep(100L, 10L),
    fun.success.runs = rep(FALSE, 10L)
  )
  expect_true(is.infinite(ERT))

  # check user-specified penalty value
  ERT = computeExpectedRunningTime(
    fun.evals = rep(100L, 10L),
    fun.success.runs = rep(FALSE, 10L),
    penalty.value = 10000L
  )
  expect_equal(ERT, 10000L)
})
