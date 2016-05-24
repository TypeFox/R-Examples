context("Testing emd_num_imfs")

test_that("bogus arguments throw error",{
  expect_error(emd_num_imfs("abc"))
  expect_error(emd_num_imfs(NULL))
  expect_error(emd_num_imfs(Inf))
  expect_error(emd_num_imfs(-1))
  expect_error(emd_num_imfs(2.1))
})
test_that("correct values are returned",{
  expect_identical(emd_num_imfs(1), 1L)
  expect_identical(emd_num_imfs(5), 2L)
  expect_identical(emd_num_imfs(10), 3L)
  expect_identical(emd_num_imfs(16), 4L)
})
