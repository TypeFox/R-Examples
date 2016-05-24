require(testthat)

context("Feasibility checks: numerical")

test_that("isFeasible",{
    expect_false(isFeasible(editmatrix(c("x +y < 0","x+y>0")), warn=FALSE))
})

test_that("isFeasible, warning",{
    expect_warning(isFeasible(editmatrix(c("x +y < 0","x+y>0")), warn=TRUE))
})

test_that("isFeasible with 0==1",{
    expect_false(isFeasible(editmatrix("0==1"), warn=FALSE))
})


context("Feasibility checks: categorical")
test_that('isFeasible with categorical data',{
    # counterintuitive example: edits are contradictory, but space of possible records is not empty.
    E <- editarray(c(
        "a %in% letters[1:3]",
        "b %in% letters[4:6]",
        "if ( a == 'a' ) b == 'd'",
        "if ( a == 'a' ) b != 'd'"
    ))
    expect_true(isFeasible(E))
    # obvious contradiction
    expect_false(isFeasible(
        editarray(c("b %in% c('x','y')","if (TRUE) b == 'x'","if(TRUE) b != 'x'"))
    ))
})







