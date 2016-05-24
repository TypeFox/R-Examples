context("Check cummean")

x <- rnorm(100)
ans <- numeric(length(x))
for (i in seq_along(x)) {
  ans[i] <- mean(x[seq_len(i)])
}

test_that("cummean returns proper answer", {
  expect_that(length(GMCM:::cummean(x)), equals(length(x)))
  expect_that(GMCM:::cummean(x[1]), equals(x[1]))
  expect_that(GMCM:::cummean(x[0]), equals(numeric(0)))
})

test_that("cummean is computed correctly", {
  expect_that(GMCM:::cummean(x), equals(ans))
})
