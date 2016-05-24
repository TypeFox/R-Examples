require(testthat)

context("FM elimination")

# test: the first example of eliminate
test_that("eliminate works fine",{

    P <- editmatrix(c(
         "4*x1 - 5*x2 - 3*x3 + z <= 0",
         "-x1 + x2 -x3 <= 2",
         "x1 + x2 + 2*x3 <= 3",
         "-x1 <= 0",
         "-x2 <= 0",
         "-x3 <= 0"))
    P1 <- eliminate(P, "x1", fancynames=TRUE)
    Ab <- matrix(c( 
        0, -0.25, -1.75, 0.25,  2,
        0, -1.25, -0.75, 0.25,  0,
        0,  2.00,  1.00, 0.00,  5,
        0,  1.00,  2.00, 0.00,  3,
        0, -1.00,  0.00, 0.00,  0,
        0,  0.00, -1.00, 0.00,  0), byrow=TRUE,nrow=6)
    H <- matrix(c(
         TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, FALSE,  TRUE, FALSE, FALSE,
        FALSE,  TRUE,  TRUE, FALSE, FALSE, FALSE,
        FALSE, FALSE,  TRUE,  TRUE, FALSE, FALSE,
        FALSE, FALSE, FALSE, FALSE,  TRUE, FALSE,
        FALSE, FALSE, FALSE, FALSE, FALSE,  TRUE), byrow=TRUE,nrow=6)
    op <- c("<=", "<=", "<=", "<=", "<=", "<=")
    expect_true(all( Ab == getAb(P1) ))
    expect_true(all( H  == getH(P1)  ))
    expect_true(all( op == getOps(P1)))
    expect_true( geth(P1) == 1 )
})







