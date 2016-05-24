context("Error message for negative weights.")

test_that("1e3 samples from cauchy",{
    expect_error(chopthin(rcauchy(1e3),N=100,eta=4,normalise=FALSE),"Negative weights not allowed")
})
