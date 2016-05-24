source('~/git/derivmkts/R/bondsimple.R')

library(testthat)

## The saved file has the parameter values used in the test along with
## the uppercased names of the functions, which have been downcased
## for this version file
#load('~/git/derivmkts/tests/testthat/option_testvalues.Rdata')
#load('option_testvalues.Rdata')
load('~/git/derivmkts/tests/testthat/option_testvalues.Rdata')
## for each function name, we will generate results believed correct
## from options.R. Then we will test against barriers.R

coupon <- 6; mat <- 20; yield <- 0.045; principal <- 100;
p2 <- 119.6451416543
p1 <- 119.5119046772
p5 <- 119.7263415797

test_that('bondpv ', {
    expect_equal(
        bondpv(coupon=coupon, mat=mat, yield=yield, principal=principal,
               freq=2),
        p2)
    expect_equal(
        bondpv(coupon=coupon, mat=mat, yield=yield, principal=principal,
               freq=1),
        p1)
    expect_equal(
        bondpv(coupon=coupon, mat=mat, yield=yield, principal=principal,
               freq=5),
        p5)
    print('bondpv okay')
}
)

test_that('bondyield ', {
    expect_equal(
        bondyield(coupon=coupon, mat=mat, price=p2,
                  principal=principal, freq=2),
        0.045, tolerance=1e-07)
    expect_equal(
        bondyield(coupon=coupon, mat=mat, price=p1,
                  principal=principal, freq=1),
        0.045, tolerance=1e-07)
    expect_equal(
        bondyield(coupon=coupon, mat=mat, price=p5,
                  principal=principal, freq=5),
        0.045, tolerance=1e-07)
    print('bondyield okay')
}
)

test_that('duration ', {
    expect_equal(
        duration(price=p2, coupon=coupon, mat=mat, 
                  principal=principal, freq=2),
        12.6353752651, tolerance=1e-07)
    expect_equal(
        duration(coupon=coupon, mat=mat, price=p1,
                  principal=principal, freq=1),
        12.8523645958, tolerance=1e-07)
    expect_equal(
        duration(coupon=coupon, mat=mat, price=p5,
                  principal=principal, freq=5),
        12.5042977451, tolerance=1e-07)
    print('duration okay')
}
)

test_that('convexity ', {
    expect_equal(
        convexity(price=p2, coupon=coupon, mat=mat, 
                  principal=principal, freq=2),
        205.9789674739, tolerance=1e-07)
    expect_equal(
        convexity(coupon=coupon, mat=mat, price=p1,
                  principal=principal, freq=1),
        207.0220342194, tolerance=1e-07)
    expect_equal(
        convexity(coupon=coupon, mat=mat, price=p5,
                  principal=principal, freq=5),
        205.3245411485, tolerance=1e-07)
    print('convexity okay')
}
)



