# tests for enigma_data fxn in taxize
context("enigma_data")

test_that("enigma_data column selection works correctly", {
  skip_on_cran()
  
  cols <- c('namelast','visitee_namelast','last_updatedby')
  res1 <- enigma_data(dataset='us.gov.whitehouse.visitor-list', select=cols)
  expect_is(res1, "enigma")
  expect_true(res1$success)
  expect_is(res1$datapath, "character")
  expect_is(res1$info, "list")
  expect_is(res1$result, "data.frame")
  expect_equal(names(res1$result), cols)
})

test_that("enigma_data works correctly for sorting data", {
  skip_on_cran()
  
  res2 <- enigma_data(dataset='us.gov.whitehouse.visitor-list', sort='+namelast')
  res2_2 <- enigma_data(dataset='us.gov.whitehouse.visitor-list', sort='-namelast')
  expect_is(res2, "enigma")
  expect_true(res2$success)
  expect_is(res2$datapath, "character")
  expect_is(res2$info, "list")
  expect_is(res2$result, "data.frame")
  expect_equal(unique(tolower(sapply(res2_2$result$namelast, function(x) substring(x, 1, 1), USE.NAMES = FALSE))), "z")
})

test_that("enigma_data works correctly to get data subset", {
  skip_on_cran()
  
  res3 <- enigma_data(dataset = 'us.gov.whitehouse.visitor-list', where = 'total_people > 5')
  expect_is(res3, "enigma")
  expect_true(res3$success)
  expect_is(res3$datapath, "character")
  expect_is(res3$info, "list")
  expect_is(res3$result, "data.frame")
  expect_more_than(min(as.numeric(res3$result$total_people)), 5)
})
