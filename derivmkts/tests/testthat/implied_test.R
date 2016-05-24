source('~/git/derivmkts/R/implied.R')

library(testthat)
## we will generally use these parameters:
#kseq <- rep(c(35, 40, 45), times=3)
#Hseq <- rep(c(35, 40, 45), each=3)
#s <- 40; v <- 0.30; r <- 0.08;  tt <- 0.25; d <- 0.02

## The saved file has the parameter values used in the test along with
## the uppercased names of the functions, which have been downcased
## for this version file
#load('~/git/derivmkts/tests/testthat/option_testvalues.Rdata')
#load('option_testvalues.Rdata')
load('~/git/derivmkts/tests/testthat/option_testvalues.Rdata')
## for each function name, we will generate results believed correct
## from options.R. Then we will test against barriers.R


for (i in imptestfns[1:2]) {
    fns <- tolower(i)
    test_that(paste(fns, 'works'), {
                  correct <- bsimpvals[, i]
                  unknown <- do.call(fns,
                                     list(s=s, k=k, r=r,
                                          tt=tt, d=d, price=prices))
                  expect_equal(correct, unknown)
              }
              )
    print(paste(i, 'okay'))
}

for (i in imptestfns[3:4]) {
    fns <- tolower(i)
    test_that(paste(fns, 'works'), {
                  correct <- bsimpvals[, i]
                  unknown <- do.call(fns,
                                     list(s=s, k=k, v=v, r=r,
                                          tt=tt, d=d, price=prices))
                  expect_equal(correct, unknown)
              }
              )
    print(paste(i, 'okay'))
}

