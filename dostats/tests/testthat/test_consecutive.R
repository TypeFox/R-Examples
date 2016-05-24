test_that("Testing seq_consecutive",{
  expect_equal(seq_consecutive(c(3,2,2,1,1,1)), c(1,2,2,3,3,3))
  expect_equal(seq_consecutive(c('a','a','a', 'b','b','b')), c(1,1,1,2,2,2))
  expect_equal(seq_consecutive('a'), 1)
  expect_error(seq_consecutive(NULL))
})