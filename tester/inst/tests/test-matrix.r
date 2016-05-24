context("Matrices")

test_that("not matrix objects return FALSE", {
  expect_that(is_matrix(1), is_false())
  expect_that(is_matrix(1:5), is_false())
  expect_that(is_matrix("string"), is_false())
  expect_that(is_matrix(iris), is_false())
  expect_that(is_matrix(iris), is_false())
  expect_that(is_matrix(list(1,2,3)), is_false())
})


test_that("not matrix objects return TRUE", {
  expect_that(is_not_matrix(1), is_true())
  expect_that(is_not_matrix(1:5), is_true())
  expect_that(is_not_matrix("string"), is_true())
  expect_that(is_not_matrix(iris), is_true())
  expect_that(is_not_matrix(list(1,2,3)), is_true())  
})


test_that("matrix return TRUE", {
  num_matrix = matrix(1:6, 2, 3)
  str_matrix = matrix(letters[1:6], 2, 3)
  expect_that(is_matrix(num_matrix), is_true())
  expect_that(is_matrix(str_matrix), is_true())
})


test_that("numeric matrix works", {
  num_matrix = matrix(1:6, 2, 3)
  str_matrix = matrix(letters[1:6], 2, 3)
  expect_that(is_numeric_matrix(num_matrix), is_true())
  expect_that(is_numeric_matrix(str_matrix), is_false())
})


test_that("string matrix works", {
  num_matrix = matrix(1:6, 2, 3)
  str_matrix = matrix(letters[1:6], 2, 3)
  expect_that(is_string_matrix(num_matrix), is_false())
  expect_that(is_string_matrix(str_matrix), is_true())
})
