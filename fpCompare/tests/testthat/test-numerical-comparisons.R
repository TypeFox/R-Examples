test_that("relational operators within tolerance", {
  x <- .5-.3
  y <- .3-.1
  expect_equal(x %==% y, TRUE)
  expect_equal(x %!=% y, FALSE)

  set.seed(123L)
  a <- jitter(1:10, 1e-3)
  b <- jitter(1:10, 1e-3)
  less.eq    <- a %<=% b
  less       <- a %<<% b
  greater.eq <- a %>=% b
  greater    <- a %>>% b
  equal      <- a %==% b
  notequal   <- a %!=% b
  expect_equal(less,     !greater)
  expect_equal(notequal, !equal)

  a <- jitter(1:10, 1e-7)
  b <- jitter(1:10, 1e-7)
  less.eq    <- a %<=% b
  less       <- a %<<% b
  greater.eq <- a %>=% b
  greater    <- a %>>% b
  equal      <- a %==% b
  notequal   <- a %!=% b
  ids <- c(5L,8L)
  expect_equal(less[ids],        !greater[ids])
  expect_equal(less[-ids],       greater[-ids])
  expect_equal(less.eq,          !greater)
  expect_equal(all(!equal[ids]), TRUE)
  expect_equal(all(equal[-ids]), TRUE)
})
