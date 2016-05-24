##
## Demo for solving an analytic centering problem with equality constraints
## (Example taken from cvxopt's userguide)
##
if(requireNamespace("numDeriv", quietly = TRUE)){
    ## Creating objective, gradient and Hessian
    f0 <- function(x) -sum(log(x))
    g0 <- function(x, func = f0) numDeriv::grad(func = func, x = x)
    h0 <- function(x, func = f0) numDeriv::hessian(func = func, x = x)
    ## equality constraint
    A <- matrix(c(1, 1, 2), nrow = 1)
    b <- matrix(1, nrow = 1)
    ## initial (feasible!) point
    x0 = c(0.25, 0.25, 0.25)
    ## solving problem
    ans <- cccp(x0 = x0, f0 = f0, g0 = g0, h0 = h0, A = A, b = b)
    ans
    getx(ans)
}
