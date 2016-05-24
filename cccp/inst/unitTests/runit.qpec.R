##
## Unit testing of Quadratic Program with equality constraints
test.QPEC <- function(){
    P <- 2 * matrix(c(2, .5, .5, 1), nrow = 2, ncol = 2)
    q <- c(1.0, 1.0)
    A <- matrix(c(1.0, 1.0), nrow = 1, ncol = 2)
    b <- 1.0
    ## Using main function of package
    ans <- cccp(P = P, q = q, A = A, b = b)
    checkTrue(ans$status == "optimal")
    checkEqualsNumeric(drop(getx(ans)), c(0.25, 0.75))
    checkEqualsNumeric(ans$state[1], 1.875)
    return()
}
