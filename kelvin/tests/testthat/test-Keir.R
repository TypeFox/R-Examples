##

context("Keir functions")

x <- seq_len(10)
ki <- Kei(x)
kr <- Ker(x)
k <- Keir(x)
k.l <- Keir(x, return.list = TRUE)

test_that("Convenience functions work properly",{
  expect_is(ki, 'matrix')
  expect_is(kr, 'matrix')
})

test_that("Primary function works properly",{
  expect_is(k, 'matrix')
  expect_is(k.l, 'list')
  expect_null(names(k))
  expect_equal(names(k.l), c('kei','ker'))
})

test_that("Primary function options work properly",{
  expect_message(Keir(x, show.scaling = TRUE))
  expect_error(Keir(c(0,x), add.tol = FALSE))
  expect_warning(k0 <- Keir(c(0,x)))
  expect_warning(k0.l <- Keir(c(0,x), return.list = TRUE))
  expect_null(names(k0))
  expect_equal(names(k0.l), c('kei','ker','zero.indices'))
})