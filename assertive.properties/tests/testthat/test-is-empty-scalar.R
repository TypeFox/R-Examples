test_that("test.is_empty.empty_list.returns_true", {
  expect_true(is_empty(list()))
})

test_that("test.is_empty.empty_vector.returns_true", {
  expect_true(is_empty(numeric()))
})

test_that("test.is_empty.non_empty_vector.returns_false", {
  expect_false(is_empty(1))
})

test_that("test.is_empty.null.returns_true", {
  expect_true(is_empty(NULL))
})

test_that("test.is_non_empty.empty_list.returns_false", {
  expect_false(is_non_empty(list()))
})

test_that("test.is_non_empty.empty_vector.returns_false", {
  expect_false(is_non_empty(numeric()))
})

test_that("test.is_non_empty.non_empty_vector.returns_true", {
  expect_true(is_non_empty(1))
})

test_that("test.is_non_empty.null.returns_false", {
  expect_false(is_non_empty(NULL))
})

test_that("test.is_non_scalar.a_scalar.returns_false", {
  expect_false(is_non_scalar(1))
})

test_that("test.is_non_scalar.a_vector.returns_true", {
  expect_true(is_non_scalar(1:2))
})

test_that("test.is_non_scalar.empty.returns_true", {
  expect_true(is_non_scalar(numeric()))
}) 

test_that("test.is_scalar.a_scalar.returns_true", {
  expect_true(is_scalar(1))
})

test_that("test.is_scalar.a_vector.returns_false", {
  expect_false(is_scalar(1:2))
})

test_that("test.is_scalar.empty.returns_false", {
  expect_false(is_scalar(numeric()))
}) 

test_that("test.is_scalar.a_single_elt_list.returns_true_with_length_metric", {
  x <- list(1:5)
  expect_true(is_scalar(x))
  expect_false(is_scalar(x, "elements"))
})


test_that("test.is_of_length.length_69.returns_true_when_length_is_69", {
  x <- 1:69
  expect_false(is_of_length(x, 68))
  expect_true(is_of_length(x, 69))
  expect_false(is_of_length(x, 70))
}) 

test_that("test.has_elements.69_elts.returns_true_when_69_elts", {
  x <- 1:69
  expect_false(has_elements(x, 68))
  expect_true(has_elements(x, 69))
  expect_false(has_elements(x, 70))
}) 

