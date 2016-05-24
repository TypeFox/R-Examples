context("listsize")

test_that("incorrect inputs raise errors", {
  bads <- list(
    letters[1:5],
    matrix(1:12, 4, 3)
  )
  
  expect_that(listsize(bads[[1]]), throws_error("A list is required"))
  expect_that(listsize(bads[[2]]), throws_error("A list is required"))
})

test_that("listsize works as expected", {
  num_list = list(1:3, 4:5, 6:9)
  str_list = list(c("a","b"), letters[1:10])
  nul_list = list(1:3, NULL, 1:10)
  
  expect_that(listsize(num_list), equals(9))
  expect_that(listsize(str_list), equals(12))
  expect_that(listsize(nul_list), equals(13))
})
