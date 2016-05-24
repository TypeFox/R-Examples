##
##  c o m p l e x s t e p . R  Complex Step Derivation
##


complexstep <- function(f, x0, h = 1e-20, ...) {
    stopifnot(is.numeric(x0), is.numeric(h), h < 1e-15)
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    try(fx0hi <- f(x0 + h*1i))
    if (inherits(fx0hi, "try-error"))
        stop("Function 'f' does not appear to accept complex arguments.")
    
    if (!is.complex(fx0hi) || !is.double(f(x0))) {
        # apply Richardson's method
        f_csd <- numderiv(f, x0, h = 0.1)$df
        warning("Some maginary part is zero: applied Richardson instead.")
    } else {
        # apply complex-step method
        f_csd <- Im(fx0hi) / h          # Im(f(x0 + h * 1i)) / h
    }
    return(f_csd)
}


grad_csd <- function(f, x0, h = 1e-20, ...) {
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    z <- f(x0)
    n <- length(x0); m <- length(z)
    if (m > 1)
        stop("Funktion 'f' does not return a scalar; call 'Jacobian_csd'.")

    G <- rep(NA, n)
    for (k in 1:n) {
        x1 <- x0
        x1[k] <- x1[k] + h * 1i
        G[k] <- Im(f(x1)) / h
    }
    return(G)
}


jacobian_csd <- function(f, x0, h = 1e-20, ...) {
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    z <- f(x0)
    n <- length(x0); m <- length(z)
    J <- matrix(NA, nrow = m, ncol = n)

    for (k in 1:n) {
        x1 <- x0
        x1[k] <- x1[k] + h * 1i
        J[ , k] <- Im(f(x1)) / h
    }
    # drop matrix dimensions if n = m = 1
    # if (m == 1 && n == 1) J <- J[,]
    return(J)
}


hessian_csd <- function(f, x0, h = 1e-20, ...) {
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    z <- f(x0)
    n <- length(x0); m <- length(z)
    if (m > 1)
        stop("Funktion 'f' does not return a scalar, as needed for Hessian.")
    
    H <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
        fi <- function(x) {
            x[i] <- x[i] + h*1i
            Im(f(x)) / h
        }
        for (j in 1:n) {
            ff <- function(x) {
                xx <- x0
                xx[j] <- x
                fi(xx)
            }
            H[i, j] <- numderiv(ff, x0[j])$df
        }
    }
    return(H)
}


laplacian_csd <- function(f, x0, h = 1e-20, ...) {
    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    z <- f(x0)
    n <- length(x0); m <- length(z)
    if (m > 1)
        stop("Funktion 'f' does not return a scalar, as needed for Lapacian.")
    
    L <- rep(NA, n)
    for (i in 1:n) {
        fi <- function(x) {
            x[i] <- x[i] + h*1i
            Im(f(x)) / h
        }
        ff <- function(x) {
            xx <- x0
            xx[i] <- x
            fi(xx)
        }
        L[i] <- numderiv(ff, x0[i])$df
    }
    return(sum(L))
}
