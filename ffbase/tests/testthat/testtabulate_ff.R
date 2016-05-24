library(testthat)

context("tabulate")

test_that("tabulate.ff works",{
   x <- c(1,1,2,3,3)
   xf <- ff(x)
   
   expect_identical(tabulate(x), tabulate.ff(xf))
   
   expect_identical(tabulate(x,2), tabulate.ff(xf,2))
   
   #expect_identical(tabulate(x,2), tabulate.ff(xf,2, TRUE)[])
})