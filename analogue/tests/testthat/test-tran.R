## Test `tran()` function

## load packages
library("testthat")
library("analogue")

context("Testing tran()")

## Imbrie and Kipp example
## load the example data
data(ImbrieKipp)

test_that("tran() handles bad inputs", {
    expect_error(tran(ImbrieKipp, method = "foo"))
})
