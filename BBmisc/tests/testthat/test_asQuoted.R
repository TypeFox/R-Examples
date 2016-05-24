context("asQuoted")

test_that("asQuoted", {
  e1 = asQuoted("x == 3")
  e2 = quote(x == 3)
  expect_equal(e1, e2)
})

