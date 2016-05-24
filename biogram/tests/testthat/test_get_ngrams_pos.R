context("Counting n-grams")

test_that("Count ngrams for different distances",{
  expect_equal(length(get_ngrams_ind(10, 2, 0)[[1]]), 9)
  expect_true(all(sapply(get_ngrams_ind(10, 3, 1), length) == 6))
})