source('~/git/derivmkts/R/barriers.R')
source('~/git/derivmkts/R/blksch.R')
#source('~/inc/R/options.R')
library(testthat)
## we will generally use these parameters:
#kseq <- rep(c(35, 40, 45), times=3)
#Hseq <- rep(c(35, 40, 45), each=3)
#s <- 40; v <- 0.30; r <- 0.08;  tt <- 0.25; d <- 0.02

## The saved file has the parameter values used in the test along with
## the uppercased names of the functions, which have been downcased
## for this version file
#load('~/git/derivmkts/tests/testthat/option_testvalues.Rdata')
load('~/git/derivmkts/tests/testthat/option_testvalues.Rdata')

## for each function name, we will generate results believed correct
## from options.R. Then we will test against barriers.R

############################################################
##  barriers
############################################################

for (i in barriertestfns) {
    testfn <- tolower(i)
    test_that(paste(testfn, 'works'), {
                  correct <- barriervals[, i]
                  unknown <- do.call(testfn,
                                     list(s=s, k=kseq, v=v, r=r,
                                          tt=tt, d=d, H=Hseq))
                  expect_equal(correct, unknown)
              }
              )
    print(paste(i, 'okay'))
}

for (i in barriertestfns2) {
    testfn <- tolower(i)
    test_that(paste(testfn, 'works'), {
                  correct <- barriervals2[, i]
                  unknown <- do.call(testfn, list(s=s, v=v, r=r,
                                             tt=tt, d=d, H=Hseq2))
                  expect_equal(correct, unknown)
              }
              )
    print(paste(testfn, 'okay'))
}

for (i in c('UR', 'DR')) {
    testfn <- tolower(i)
    test_that(paste(testfn, 'perpetual works'), {
        correct <- barriervals3[, i]
##        print(correct)
                  unknown <- do.call(testfn,
                                     list(s=s, v=v, r=r, tt=tt,
                                          d=d, H=Hseq2, perpetual=TRUE)
                                     )
##        print(unknown)
##                  print(paste(s, v, r, tt, d, H))
                  expect_equal(correct, unknown)
              }
              )
    print(paste(testfn, 'okay'))
}
    
