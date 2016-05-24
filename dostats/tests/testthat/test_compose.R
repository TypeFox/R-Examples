test_that("Testing compose", {
  expect_is(compose(any, is.na), 'function')
  expect_error(compose(any, isna))
  expect_true(compose(any,is.na)(NA))
  expect_true((any%.%is.na)(NA))
})
test_that("redirf", {
exf <- function(..., sep=" ", collapse=NULL)
{
    paste(..., sep=sep, collapse=collapse)
}
testf <- redirf(paste)
identical(deparse(exf), deparse(testf), F,F,F,F)
})


