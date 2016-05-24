context('selected_columns')

test_that('selected_columns returns the correct result', {
  expect_equal(autovarCore:::selected_columns(7), c(1, 2, 3))
  expect_equal(autovarCore:::selected_columns(15), c(1, 2, 3, 4))
  expect_equal(autovarCore:::selected_columns(13), c(1, 3, 4))
  expect_equal(autovarCore:::selected_columns(511), 1:9)
})
