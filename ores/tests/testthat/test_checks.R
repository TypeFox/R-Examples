context("Test that check_* functions work")

test_that("Checking for damaging content works",{
  
  result <- check_damaging("enwiki", 34854345)
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_equal(complete.cases(result), TRUE)
})

test_that("Checking for good-faith content works",{
  
  result <- check_goodfaith("enwiki", 34854345)
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_equal(complete.cases(result), TRUE)
})

test_that("Checking for reverted content works",{
  
  result <- check_reverted("enwiki", 34854345)
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_equal(complete.cases(result), TRUE)
})

test_that("Error handling works",{
  
  result <- check_goodfaith("enwiki", 9999999999)
  expect_equal(class(result), "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 5)
  expect_equal(complete.cases(result), FALSE)
})