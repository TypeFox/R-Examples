context("anylength")
test_that("vector has correct length", {
  a <- c(1,2,3)
  expect_equal(anylength(a), 3)
})

test_that("matrix has correct length", {
  m <- matrix(1:10, ncol=2)
  expect_equal(anylength(m), 5)
})

test_that("data.frame has correct length", {
  df <- data.frame(x1=1:10, x2=1:10)
  expect_equal(anylength(df), 10)
})


context("anynames")
test_that("named matrix has correct names", {
  m <- matrix(c(1,2,3,4,5,6), ncol=2)
  anynames(m) <- c('d','e')
  expect_equal(anynames(m), c('d', 'e'))
})

test_that("named vector has correct names", {
  v <- c(a=1,b=2,c=3,d=4,e=5)
  expect_equal(anynames(v), c('a', 'b', 'c', 'd', 'e'))
})

test_that("a named list has correct names", {
  l <- list(a=1,b=2,c=3,d=4,e=5)
  expect_equal(anynames(l), c('a', 'b', 'c', 'd', 'e'))
})

test_that("a named data.frame has correct names", {
  df <- data.frame(a=1:10, b=1:10,c=1:10,d=1:10,e=1:10)
  expect_equal(anynames(df), c('a', 'b', 'c', 'd', 'e'))
})


context("anytypes")
test_that("A named data.frame has the correct types", {
  a <- data.frame(a=c(1,2,3), b=c("larry","mo","curly"), c=c(TRUE,FALSE,TRUE))
  ts <- anytypes(a)
  expect_equal(names(ts), c("a","b","c"))
  names(ts) <- NULL
  expect_equal(ts, c('numeric','factor','logical'))
})

test_that("An unnamed data.frame has the correct types", {
  a <- data.frame(c(1,2,3), c("larry","mo","curly"), c(TRUE,FALSE,TRUE))
  ts <- anytypes(a)
  # The data.frame will fill this in
  expect_true(! is.null(names(ts)))
  names(ts) <- NULL
  expect_equal(ts, c('numeric','factor','logical'))
})


context("is.bad")
test_that("A list is handled correctly", {
  a <- list(a=1:3, b=NULL, c=NA, d='foo')
  e <- list(a=rep(FALSE,3), b=TRUE, c=TRUE, d=FALSE)
  expect_equal(is.bad(a), e)
})

test_that("A vector with NAs is handled correctly", {
  a <- c(1,NA,3)
  expect_equal(is.bad(a), c(FALSE,TRUE,FALSE))
})

test_that("A data.frame with NAs is handled correctly", {
  a <- data.frame(a=1:3, b=NA)
  e <- matrix(c(rep(FALSE,3), rep(TRUE,3)), ncol=2)
  colnames(e) <- c('a','b')
  expect_equal(is.bad(a), e)
})

test_that("A data.frame that is empty is handled correctly", {
  a <- data.frame(a=NULL, b=NULL)
  expect_true(is.bad(a))
})

test_that("A matrix with NAs is handled correctly", {
  a <- matrix(c(1:3, NA), ncol=2)
  e <- matrix(c(rep(FALSE,3), TRUE), ncol=2)
  expect_equal(is.bad(a), e)
})


context("is.empty")
test_that("A non-empty vector resolves to FALSE", {
  a <- c(1,2,3)
  expect_false(is.empty(a))
})

test_that("An empty vector resolves to TRUE", {
  a <- c()
  expect_true(is.empty(a))
})

test_that("A non-empty list resolves to FALSE", {
  a <- list(a=1,2,3)
  expect_false(is.empty(a))
})

test_that("An empty list resolves to TRUE", {
  a <- list()
  expect_true(is.empty(a))
})

test_that("A non-empty data.frame resolves to FALSE", {
  a <- data.frame(a=1:3,b=2,c=3)
  expect_false(is.empty(a))
})

test_that("An empty data.frame resolves to TRUE", {
  a <- data.frame(a=NULL, b=NULL)
  expect_true(is.empty(a))
})
