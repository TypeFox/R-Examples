test_that("test.is_in_future.then_now_soon.returns_true_in_past", {
  x <- Sys.time() + c(-1, 5, NA)
  expected <- c(FALSE, TRUE, NA)
  expect_equal(strip_attributes(actual <- is_in_future(x)), expected)
  expect_named(actual)
  expect_equal(cause(actual), noquote(c("in past", "", "missing")))
})

test_that("test.is_in_past.then_now_soon.returns_true_in_past", {
  x <- Sys.time() + c(-1, 5, NA)
  expected <- c(TRUE, FALSE, NA)
  expect_equal(strip_attributes(actual <- is_in_past(x)), expected)
  expect_named(actual)
  expect_equal(cause(actual), noquote(c("", "in future", "missing")))
}) 
