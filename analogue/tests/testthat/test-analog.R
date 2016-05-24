## Test `analog()` function

## load packages
library("testthat")
library("analogue")

context("Testing analog()")


## Imbrie and Kipp example
## load the example data
data(ImbrieKipp)
data(SumSST)
data(V12.122)

## merge training and test set on columns
dat <- join(ImbrieKipp, V12.122, verbose = TRUE)

## extract the merged data sets and convert to proportions
ImbrieKipp <- tran(dat[[1]], method = "pcent2prop")
V12.122 <- tran(dat[[2]] / 100, method = "pcent2prop")

## Imbrie and Kipp foraminfera sea-surface temperature

test_that("analog() works & returns correct object", {
    ## analog matching between SWAP and RLGH core
    ik.analog <- analog(ImbrieKipp, V12.122, method = "chord")
    expect_is(ik.analog, "analog")
})
