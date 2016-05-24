##

context("Beir functions")

x <- seq_len(10)
bi <- Bei(x)
br <- Ber(x)
b <- Beir(x)
b.l <- Beir(x, return.list = TRUE)
b2 <- Beir(x, nSeq=2)

test_that("Convenience functions work properly",{
  expect_is(bi, 'matrix')
  expect_is(br, 'matrix')
})

test_that("Primary function works properly",{
  expect_is(b, 'matrix')
  expect_is(b.l, 'list')
  expect_null(names(b))
  expect_equal(names(b.l), c('bei','ber'))
})
