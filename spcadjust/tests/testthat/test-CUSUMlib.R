context("Markov Approximation for properties of CUSUM charts")


test_that("Extreme case: almost constant updates",{
    skip("Exreme case that is not expected to work")
    expect_equal(ARL_CUSUM_Markovapprox(3,pobs=function(x)pnorm(x,mean=0,sd=0.0025)),
                 ARL_CUSUM_Markovapprox(3,gridpoints=500,pobs=function(x)pnorm(x,mean=0,sd=0.0025)),tolerance=0.1)
})


test_that("Hitting probability calibration - skew distr with large positive components",{
    skip("Exreme case that is not expected to work")
    set.seed(347892)
    X <-  rexp(100)
    eobs <- ecdf((X-mean(X)-0.05)/sd(X))
    expect_equal(calibratehitprob_Markovapprox(n=10000,pobs=eobs),
                 calibratehitprob_Markovapprox(n=10000,pobs=eobs,gridpoints=200),tolerance=0.1)
})
