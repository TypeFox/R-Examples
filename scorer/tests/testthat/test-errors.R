library(scorer)
context("errors")

n <- 10000
x <- 1:n
test_that("e() produces correct output.", {
  expect_equal(e(x, x), rep(0, n))
  expect_equal(e(0, 1), -1)
  expect_equal(e(1, 0), 1)
  expect_equal(e(FALSE, TRUE), -1)
  expect_equal(e(TRUE, FALSE), 1)
})

test_that("e() produces correct output classes and types.", {
  expect_is(e(runif(n), runif(n)), "numeric")
})

test_that("e() raises correct messages, warnings, and errors.", {
  expect_warning(e(x, x))
  expect_error(e("a", "b"))
})

test_that("ae() produces correct output.", {
  expect_equal(ae(x, x), rep(0, n))
  expect_equal(ae(0, 1), 1)
  expect_equal(ae(1, 0), 1)
  expect_equal(ae(FALSE, TRUE), 1)
  expect_equal(ae(TRUE, FALSE), 1)
})

test_that("ae() produces correct output classes and types.", {
  expect_is(ae(runif(n), runif(n)), "numeric")
})

test_that("ae() raises correct messages, warnings, and errors.", {
  expect_warning(ae(x, x))
  expect_error(ae("a", "b"))
})

test_that("ape() produces correct output.", {
  expect_equal(ape(x, x), rep(0, n))
  expect_equal(ape(0, 1), Inf)
  expect_equal(ape(1, 0), 1)
  expect_equal(ape(FALSE, TRUE), Inf)
  expect_equal(ape(TRUE, FALSE), 1)
})

test_that("ape() produces correct output classes and types.", {
  expect_is(ape(runif(n), runif(n)), "numeric")
})

test_that("ape() raises correct messages, warnings, and errors.", {
  expect_warning(ape(x, x))
  expect_error(ape("a", "b"))
})

test_that("mape() produces correct output.", {
  expect_equal(mape(x, x), 0)
  expect_equal(mape(0, 1), Inf)
  expect_equal(mape(1, 0), 1)
  expect_equal(mape(FALSE, TRUE), Inf)
  expect_equal(mape(TRUE, FALSE), 1)
})

test_that("mape() produces correct output classes and types.", {
  expect_is(mape(runif(n), runif(n)), "numeric")
})

test_that("mape() raises correct messages, warnings, and errors.", {
  expect_warning(mape(x, x))
  expect_error(mape("a", "b"))
})

test_that("pe() produces correct output.", {
  expect_equal(pe(x, x), rep(0, n))
  expect_equal(pe(0, 1), -Inf)
  expect_equal(pe(1, 0), 1)
  expect_equal(pe(FALSE, TRUE), -Inf)
  expect_equal(pe(TRUE, FALSE), 1)
})

test_that("pe() produces correct output classes and types.", {
  expect_is(pe(runif(n), runif(n)), "numeric")
})

test_that("pe() raises correct messages, warnings, and errors.", {
  expect_warning(pe(x, x))
  expect_error(pe("a", "b"))
})


test_that("rmse() produces correct output.", {
  expect_equal(rmse(x, x), 0)
  expect_equal(rmse(0, 1), 1)
  expect_equal(rmse(1, 0), 1)
  expect_equal(rmse(2, 0), 2)
  expect_equal(rmse(0, 2), 2)
  expect_equal(rmse(FALSE, TRUE), 1)
  expect_equal(rmse(TRUE, FALSE), 1)
})

test_that("rmse() produces correct output classes and types.", {
  expect_is(rmse(runif(n), runif(n)), "numeric")
})

test_that("rmse() raises correct messages, warnings, and errors.", {
  expect_warning(rmse(x, x))
  expect_error(rmse("a", "b"))
})
