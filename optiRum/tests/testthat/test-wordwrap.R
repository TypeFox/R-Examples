context("wordwrap")
library(plyr)
test_that("wordwrap format always adds \n", {
  expect_equal(wordwrap("a balloon"), "a\nballoon")
  expect_equal(wordwrap("a"), "a")
})
