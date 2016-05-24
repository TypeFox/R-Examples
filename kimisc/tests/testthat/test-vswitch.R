context("vswitch")

test_that("character", {
  expect_error(vswitch(LETTERS, 1, 2, 3), "unnamed")
  expect_equal(vswitch(LETTERS), rep(NA, length(LETTERS)))
  expect_equal(vswitch(LETTERS, A=1, B=2, 3),
               c(1:3, rep(3, length(LETTERS) - 3)))
  expect_equal(vswitch(LETTERS, A=1, B=1 + 1, 3, D=2 * 2),
               c(1:4, rep(3, length(LETTERS) - 4)))
  expect_equal(vswitch(LETTERS, ZZ=stop()), rep(NA, length(LETTERS)))
  expect_error(vswitch(LETTERS, A=stop("Found A")), "Found")
})

test_that("default vector", {
  expect_equal(vswitch(LETTERS, A=1, B=2, letters),
               c(1:2, rep(list(letters), 24)))
})

test_that("numeric", {
  skip("NYI")
})

test_that("factor", {
  expect_error(vswitch(factor(LETTERS)), "EXPR")
})

test_that("logical", {
  expect_error(vswitch(TRUE, "EXPR"))
})
