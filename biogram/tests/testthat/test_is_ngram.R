context("n-gram validation")

test_that("Real n-grams",{
  
  #n-gram with position information
  expect_true(is_ngram("1_1.1.1_0.0"))
  
  #n-gram without position information
  expect_true(is_ngram("1.1.1_0.0"))
  
  #False n-grams
  
  #n-gram with position information but no distance
  expect_false(is_ngram("1_1.1.1"))
})
