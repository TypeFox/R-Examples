##
##  r o m b e r g . R  Romberg Integration
##


romberg <- function(f, a, b, maxit = 25, tol = 1e-12, ...)
{
    stopifnot(is.numeric(a), is.numeric(b), length(a) == 1, length(b) == 1)
    tol <- abs(tol)
    if (a == b) return(list(value = 0, iter = 0, rel.error = 0))
    if (a > b)  return(-1 * romberg(f, b, a, tol = tol, ...))

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    eps <- .Machine$double.eps
    if (!is.finite(f(a))) a <- a + eps * sign(b-a)
    if (!is.finite(f(b))) b <- b - eps * sign(b-a)

    I <- matrix(0, nrow = maxit+1, ncol = maxit+1)
    n <- 1; iter <- 0; err <- 1.0
    while (err > tol && iter < maxit) {
        iter <- iter+1

        # inner trapezoid rule with correction term
        n <- 2 * n
        h <- (b-a) / n
        S <- f(a)
        for (i in 1:(n-1)) {
            xi <- a + h * i
            S <- S + 2*f(xi)
        }
        S <- (S + f(b)) * h/2
        # S <- S - ((f(b)-f(xi)) - (f(a+h)-f(a))) * h/12.0  # truncation errors
        I[iter+1, 1] <- S

        # Richardson approximation
        for (k in 2:(iter+1)) {
            j <- 2+iter-k
            I[j,k] <- (4^(k-1)*I[j+1,k-1] - I[j,k-1]) / (4^(k-1)-1)
        }
        err <- abs((I[1,iter+1] - I[2,iter]) / I[1,iter+1])
    }

    if (iter == maxit)
        warning("Maximum number of iterations has been reached.")
    if (err < tol) err <- tol 

    return(list(value = I[1, iter+1], iter = iter, rel.error = err))
}
