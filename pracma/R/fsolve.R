##
##  f s o l v e . R
##


fsolve <- function(f, x0, J = NULL,
                   maxiter = 100, tol = .Machine$double.eps^(0.5), ...) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.")
    x0 <- c(x0)

    # Prepare objective function and its Jacobian
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    n <- length(x0)
    m <- length(f(x0))

    if (!is.null(J)) {
        Jun <- match.fun(J)
        J <- function(x) J(x, ...)
    } else {
        J <- function(x) jacobian(f, x)
    }

    if (m == n) {
        sol = broyden(f, x0, J0 = J(x0), maxiter = maxiter, tol = tol)
        xs <- sol$zero; fs <- f(xs)
    } else {
        sol <- gaussNewton(x0, f, Jfun = J, maxiter = maxiter, tol = tol)
        xs <- sol$xs; fs <- sol$fs
        if (fs > tol)
            warning("Minimum appears not to be a zero -- change starting point.")
    }

    return(list(x = xs, fval = fs))
}


fzsolve <- function(fz, z0) {
    if (length(z0) == 0) return(c())
    if (length(z0) > 1) {
        warning("Argument 'z0' has length > 1, first component taken.")
        z0 <- z0[1]
    }

    if (is.numeric(z0)) {
        x0 <- c(z0, 0)
    } else if (is.complex(z0)) {
        x0 <- c(Re(z0), Im(z0))
    } else
        stop("Argument 'z0' must be a real or complex number.")

    fz <- match.fun(fz)
    fn <- function(x) {
        z <- x[1] + x[2]*1i
        Z <- fz(z)
        if (is.complex(Z))      Z <- c(Re(Z), Im(Z))
        else if (is.numeric(Z)) Z <- c(Z, 0)
        else                    Z <- c(NA, NA)
        return(Z)
    }
    x <- broyden(fn, x0)$zero
    return(x[1] + x[2]*1i)
}
