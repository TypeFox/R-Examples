library(profilr)
context("concat_rank_value_occurrence()")

test_that("concat_rank_value_occurrence of character concatenates properly ", {
  expect_equal(concat_rank_value_occurrence("a", "b", "c"), "(a) b [c]\n")
  expect_equal(concat_rank_value_occurrence(1, 2, 3), "(1) 2 [3]\n")
  expect_equal(concat_rank_value_occurrence(5, TRUE, 10), "(5) TRUE [10]\n")
})
