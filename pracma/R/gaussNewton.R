##
##  g a u s s N e w t o n . R  Gauss-Newton Function Minimization
##


gaussNewton <- function(x0, Ffun, Jfun = NULL,
                         maxiter =100, tol = .Machine$double.eps^(1/2), ...) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.")
    x0 <- c(x0)
    n <- length(x0)
    if (n == 1)
        stop("Function is univariate -- use a different approach.")

    Fun <- match.fun(Ffun)
    F <- function(x) Fun(x, ...)
    m <- length(F(x0))

    # Define J as Jacobian of F, and f as sum(F_i)
    f <- function(x) sum(F(x)^2)
    if (is.null(Jfun)) {
        J <- function(x) jacobian(F, x)
    } else {
        Jun <- match.fun(Jfun)
        J <- function(x) Jun(x, ...)
    }
    # The gradient g of f (!) computed through the Jacobian
    g <- function(x) 2 * t(J(x)) %*% F(x)

    # First iteration
    xk <- x0
    fk <- f(xk)
    Jk <- J(xk)
    gk <- g(xk)  # 2 * t(Jk) %*% F(xk)
    Hk <- 2 * t(Jk) %*% Jk + tol * diag(1, n)

    # Line search step
    dk <- -inv(Hk) %*% gk
    if (!all(is.finite(dk))) dk <- .mdhess(Hk, gk)
    ak <- softline(xk, dk, f, g)
    adk <- ak * dk
    err <- Norm(adk)
    
    k <- 1
    while (err > tol && k < maxiter) {
        xk <- xk + adk
        fl <- f(xk)
        Jk <- J(xk)
        gk <- g(xk)
        Hk <- 2 * t(Jk) %*% Jk + tol * diag(1, n)

        dk <- -inv(Hk) %*% gk
        if (!all(is.finite(dk))) dk <- .mdhess(Hk, gk)
        ak <- softline(xk, dk, f, g)
        adk <- ak * dk
        err <- abs(fl - fk)

        # Prepare next iteration
        fk <- fl
        k <- k + 1
    }

    if (k >= maxiter)
        warning("Maximum number of iterations reached -- no convergence.")
    xs <- c(xk + adk)
    fs <- f(xk + adk)
        
    return(list(xs = xs, fs = fs, niter = k, relerr = err))
}


.mdhess <- function(H, gk) {
    n <- length(gk)

    # Matthews-Davies algorithm
    L <- D <- matrix(0, n, n)
    h00 <- if (H[1, 1] > 0) H[1, 1] else 1
    for (k in 2:n) {
        m <- k-1
        L[m, m] <- 1
        if (H[m, m] <= 0) H[m, m] <- h00
        for (i in k:n) {
            L[i, m] <- -H[i, m] / H[m, m]
            H[i, m] <- 0
            for (j in k:n) {
                H[i, j] <- H[i, j] + L[i, m] * H[m, j]
            }
        }
        if (H[k, k] > 0 && H[k, k] < h00) h00 <- H[k, k]
    }
    L[n, n] <- 1
    if (H[n,n] <= 0) H[n, n] <- h00
    for (i in 1:n)
        D[i, i] <- H[i, i]

    # Determine direction vector
    yk <- -L %*% gk
    dk <- t(L) %*% diag(1/diag(D)) %*% yk

    return(dk)
}
