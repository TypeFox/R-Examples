context("Feature Class")

test_that("Creating features works", {
  feat <- feature_class$new("abc", 5)
  expect_true(is.feature(feat))
  expect_false(is.feature("blub"))
  expect_equal(feat$get_parameters(), list())
})


test_that("adding parameter works", {
  feat <- R6::R6Class("Feature_test", inherit = feature_class,
    public = list(
      initialize = function() {}, #nolint
      add_parameter = function(...) private$add_parameter(...)
    )
  )$new()
  expect_equal(feat$add_parameter(5), "par(5)")
  expect_equal(feat$add_parameter("17"), "par(17)")
  expect_equal(feat$add_parameter(par_expr(tau)), "par(tau)")
  expect_equal(length(feat$get_parameters()), 1)
  expect_equal(feat$add_parameter(par_const(5)), "par(5)")

  expect_error(feat$add_parameter(1:2))
  expect_error(feat$add_parameter(NA))
  expect_error(feat$add_parameter(NULL))

  expect_equal(feat$add_parameter(NA, FALSE), NA)
  expect_equal(feat$add_parameter(NULL, FALSE), NA)
})


test_that("print parameter works", {
  expect_equal(print_par("par(5)"), "`5`")
  expect_equal(print_par("par(2 * theta)"), "`2 * theta`")
})
