context("Testing distname-utils \n")

test_that("distnames are in alphabetical order", {
  
  expect_true(!is.unsorted(get_distnames()))
})