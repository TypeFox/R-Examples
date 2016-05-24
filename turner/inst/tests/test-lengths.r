context("lengths")

test_that("incorrect inputs raise errors", {
  bads <- list(
    letters[1:5],
    matrix(1:12, 4, 3)
  )
  
  expect_that(lengths(bads[[1]]), throws_error("A list is required"))
  expect_that(lengths(bads[[2]]), throws_error("A list is required"))
})

test_that("lengths works as expected", {
  num_list = list(1:3, 4:5, 6:9)
  str_list = list(c("a","b"), letters[1:10])
  nul_list = list(1:3, NULL, 1:10)
  
  expect_that(lengths(num_list), equals(c(3, 2, 4)))
  expect_that(lengths(num_list, 'list'), equals(list(3, 2, 4)))
  expect_that(lengths(str_list), equals(c(2, 10)))
  expect_that(lengths(str_list, 'list'), equals(list(2, 10)))  
  expect_that(lengths(nul_list), equals(c(3, 0, 10)))
  expect_that(lengths(nul_list, 'list'), equals(list(3, 0, 10)))  
})
