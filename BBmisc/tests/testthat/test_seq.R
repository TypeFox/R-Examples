context("seq")

test_that("seq", {
  expect_equal(seq_row(iris), seq_len(nrow(iris)))
  expect_equal(seq_col(iris), seq_len(ncol(iris)))
})
