# :vim set filetype=R
context("range_for")
test_that("range_for works for a sequence with repeating values", {
  sequence <- c(seq(1, 25), c(25, 25), seq(26, 50))
  y <- data.frame(min=25, max=27)
  expect_equal(range_for(25, sequence), y)
})

test_that("Series can not be a 2-d array", {
  x <- matrix(1:10, ncol=2)
  expect_error(range_for(2, x), "No valid function for")
})


context("samplerange")
test_that("Window can not be greater than the length of x", {
  x <- rnorm(10)
  expect_error(samplerange(x, 20, 50), "No valid function for") 
})
