context("Parameter Ranges")

test_that("Ranges can be initialized", {
  par_ranges_class$new(matrix(1:4, 2, 2))
  par_ranges_class$new(matrix(1:6, 3, 2))
  par_ranges_class$new(matrix(-3:2, 3, 2))
  expect_error(par_ranges_class$new(1:6))
  expect_error(par_ranges_class$new(matrix(1:6, 2, 3)))
  expect_error(par_ranges_class$new(matrix(6:1, 3, 2)))
})


test_that("Parameter normalization works", {
  par_range <- par_ranges_class$new(matrix(-3:2, 3, 2))
  expect_true(all(par_range$normalize(-3:-1) == c(0, 0, 0)))
  expect_true(all(par_range$normalize(0:2) == c(1, 1, 1)))
  expect_true(all(par_range$normalize(-2:0) > 0))
  expect_true(all(par_range$normalize(-2:0) < 1))
})


test_that("Parameter denormalization works", {
  par_range <- par_ranges_class$new(matrix(-3:2, 3, 2))
  expect_equivalent(par_range$denormalize(c(0, 0, 0)), -3:-1)
  expect_equivalent(par_range$denormalize(c(1, 1, 1)), 0:2)  
  expect_true(all(par_range$normalize(rep(.5, 3)) > -3))
  expect_true(all(par_range$normalize(rep(.5, 3)) < 2))  
})
