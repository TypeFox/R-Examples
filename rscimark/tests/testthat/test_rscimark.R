test_that("rscimark", {
  expect_scimark_output = function(x) {
    expect_true(is.numeric(x))
    expect_true(length(x) == 6L)
    expect_true(!any(is.na(x) | is.nan(x)))
    expect_true(!is.null(names(x)))
  }

  y = rscimark(minimum.time = 0.2)
  expect_scimark_output(y)

  y = rscimark(large = TRUE, minimum.time = 0.2)
  expect_scimark_output(y)
})
