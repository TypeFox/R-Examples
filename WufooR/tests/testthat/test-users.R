library(testthat)
library(WufooR)

options(Wufoo_Name = "johnmalc", Wufoo_API = "F1QH-Q64B-BSBI-JASJ")

context("Users")

test_that("User request returns 17 rows, always", {
  userDB <- user_info()
  expect_equal(dim(userDB)[2], 17)
})
