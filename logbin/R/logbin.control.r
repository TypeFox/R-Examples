logbin.control <- function (bound.tol = 1e-06, epsilon = 1e-08, maxit = 10000, trace = 0) 
{
    if (!is.numeric(bound.tol) || bound.tol <= 0)
        stop("value of 'bound.tol' must be > 0")
    if (!is.numeric(epsilon) || epsilon <= 0) 
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0) 
        stop("maximum number of iterations must be > 0")
    list(bound.tol = bound.tol, epsilon = epsilon, maxit = maxit, trace = trace)
}