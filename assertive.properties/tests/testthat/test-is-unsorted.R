test_that("test.is_unsorted.a_sorted_vector.returns_false", {
  expect_false(is_unsorted(1:3))
})

test_that("test.is_unsorted.an_unsorted_vector.returns_true", {
  expect_true(is_unsorted(c(1, 3, 2)))
})

test_that("test.is_unsorted.an_weakly_unsorted_vector.returns_strictly", {
  expect_false(is_unsorted(c(1, 1, 2)))
  expect_true(is_unsorted(c(1, 1, 2), strictly = TRUE))
})
