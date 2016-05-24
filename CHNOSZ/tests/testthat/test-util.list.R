context("util.list")

test_that("which.pmax() properly applies attributes, and also works for lists of length 1", {
  testlist <- list(a=matrix(c(1,2,3,4)), b=matrix(c(4,3,2,1)))
  testattr <- attributes(testlist[[1]])
  expect_equal(attributes(which.pmax(testlist)), testattr)
  expect_equal(as.numeric(which.pmax(testlist[1])), c(1, 1, 1, 1))
})
