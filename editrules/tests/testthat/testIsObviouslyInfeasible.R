library(testthat)

context("Obvious infeasibility")

test_that("Obvious infeasibility is detected",{
    expect_true(isObviouslyInfeasible(editmatrix("0*x == 1")))    
    expect_true(isObviouslyInfeasible(editmatrix("0*x < -1")))
    expect_true(isObviouslyInfeasible(editmatrix("1e-12*x <= -1")))    
    expect_true( isObviouslyInfeasible(editmatrix("0*x  < 0")))
    expect_true( isObviouslyInfeasible(editmatrix("0*x  < 1e-12")))
    expect_true( isObviouslyInfeasible(editmatrix("1e-12*x  < 1e-12")))
    expect_false(isObviouslyInfeasible(editmatrix("0*x  <= 0")))
    expect_false(isObviouslyInfeasible(editmatrix("x  <= 0")))
    expect_false(isObviouslyInfeasible(editmatrix("0*x  <=1e-12")))
})




