context("itostr")

test_that("itostr", {
  x = 0:200
  base = c(2, 10, 16, 24, 32, 36)

  for (b in base) {
    res = itostr(x, b)
    expect_true(is.character(res))
    expect_true(length(x) == length(res))
    expect_equal(strtoi(res, b), x)
  }

  expect_error(itostr(1, 0))
  expect_error(itostr(1, 1))
  expect_error(itostr(-1, 2))
})
