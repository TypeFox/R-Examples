broyden <- function(Ffun, x0, J0 = NULL, ...,
                    maxiter = 100, tol = .Machine$double.eps^(1/2)) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric (row or column) vector.")
    fun <- match.fun(Ffun)
    F <- function(x) fun(x, ...)
    y0 <- F(x0)
    if (length(x0) != length(y0))
        stop("Function 'F' must be 'square', i.e. from R^n to R^n .")

    # Compute once the Jacobian and its inverse
    if (is.null(J0)) {
        A0 <- jacobian(F, x0)
    } else {
        A0 <- J0
    }

    B0 <- inv(A0)
    if (any(is.infinite(B0)))
        B0 <- diag(length(x0))

    # Secant-like step in Broyden's method
    xnew <- x0 - B0 %*% y0
    ynew <- F(xnew)

    k <- 1
    while (k < maxiter) {
        s <- xnew - x0
        d <- ynew - y0
        if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < tol)  break

        # Sherman-Morrison formula
        B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0 / c(t(s) %*% B0 %*% d)

        # Update for next iteration
        x0 <- xnew
        xnew <- xnew - B0 %*% ynew
        y0 <- ynew
        ynew <- F(xnew)

        k <- k + 1
    }

    if (k >= maxiter)
        warning(paste("Not converged: Max number of iterations reached."))

    fnew <- sqrt(sum(ynew^2))
    return(list(zero = c(xnew), fnorm = fnew, niter = k))
}
