library(testthat)
library(quhomology)
context("Homology")

test_that("homology calculates correctly:",{
  expect_equal(homology(2,3,T,T),c(1,1,1,1))
  expect_equal(homology(2,3,F,T),c(1,1,1,1,1,1,0))
  expect_equal(degenerate_homology(3,3,T),c(1,1,1,1,1,1,1,1,1,1,1,1,0))
})


context("S Test")
result <- rep(T,4)
names(result) <- c("S permutation","f permutation", "g permutation", "Yang-Baxter")
test_that("S Test calculates correctly:",{
  expect_equal(quhomology::S_test(3,T),result)
})
