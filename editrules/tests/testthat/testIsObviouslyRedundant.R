
require(testthat)
context("Obvious redundancy")

test_that("Obviously redundant rows are detected",{
    expect_true(isObviouslyRedundant(editmatrix("0*x == 1e-12")))
    expect_true(isObviouslyRedundant(editmatrix("0*x <= 0")))
    expect_true(isObviouslyRedundant(editmatrix("0*x < 1")))
    expect_false(isObviouslyRedundant(editmatrix("0*x < 0")))
    expect_true(isObviouslyRedundant(editmatrix("0*x <= 1")))
    expect_true(isObviouslyRedundant(editmatrix("1e-12*x < 1")))
    
    expect_equal(isObviouslyRedundant(editmatrix(c("0*x <= 0"
                                                  , "y < 0"
                                                  )
                                                )
                                     )
                , c(TRUE,FALSE)
                )
})




