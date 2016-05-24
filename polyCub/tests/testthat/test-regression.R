context("Regression tests")

test_that("isotropic cubature can handle control list for integrate()", {
    data("letterR", package="spatstat", envir=environment())
    f <- function (s) (rowSums(s^2)+1)^-2
    ## previosly, passing control arguments did not work
    int1 <- polyCub.iso(letterR, f, center=c(0,0), control=list(rel.tol=1e-3))
    int2 <- polyCub.iso(letterR, f, center=c(0,0), control=list(rel.tol=1e-8))
    ## results are almost identical
    expect_that(int1, equals(int2, tolerance=1e-3))
    expect_that(int1, not(is_identical_to(int2)))
})
