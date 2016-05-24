library(testthat)
library(ScrabbleScore)

#no sense in testing every word, just setting up framework for future bugs
test_that("is.twl06.word return TRUE for words in the list and FALSE for those not",{
  expect_true(is.twl06.word("zzz"))
  expect_false(is.twl06.word("zzzz"))
  
})