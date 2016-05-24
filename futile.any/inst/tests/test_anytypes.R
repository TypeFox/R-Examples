context("anytypes")

test_that("unnamed vector has correct type", {
  a <- c(1,2,3)
  expect_that(anytypes(a), equals("numeric"))
})

test_that("named vector has correct type", {
  a <- c(1,2,3)
  names(a) <- c('a','b','c')
  expect_that(anytypes(a), equals("numeric"))
})

test_that("named data.frame has correct type", {
  a <- data.frame(a=c(1,2,3), b=c("larry","mo","curly"), c=c(TRUE,FALSE,TRUE))
  ts <- anytypes(a)
  expect_that(names(ts), equals(c("a","b","c")))
  names(ts) <- NULL
  expect_that(ts, equals(c('numeric','factor','logical')))
})

test_that("unnamed data.frame has correct type", {
  a <- data.frame(c(1,2,3), c("larry","mo","curly"), c(TRUE,FALSE,TRUE))
  ts <- anytypes(a)
  # The data.frame will fill this in
  expect_that(! is.null(names(ts)), is_true())
  names(ts) <- NULL
  expect_that(ts, equals(c('numeric','factor','logical')))
})

