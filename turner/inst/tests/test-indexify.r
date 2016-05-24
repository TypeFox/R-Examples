context("indexify lists")

test_that("incorrect inputs raise errors", {
  bads <- list(
    c(1:4, 4:7),
    matrix(1:12, 4, 3)
  )
  
  expect_that(indexify(bads[[1]]), throws_error())
  expect_that(indexify(bads[[2]]), throws_error())
})

test_that("indexify works as expected", {
  num_list = list(rnorm(3), rnorm(5))
  str_list = list(c("a","b","c"), c("d", "e"), c("f"))
  
  expect_that(indexify(num_list), equals(c(rep(1,3), rep(2,5))))
  expect_that(indexify(num_list, 'list'), equals(list(rep(1,3), rep(2,5))))
  expect_that(indexify(str_list), equals(c(1,1,1,2,2,3)))
  expect_that(indexify(str_list, 'list'), equals(list(c(1,1,1), c(2,2), 3)))
})
