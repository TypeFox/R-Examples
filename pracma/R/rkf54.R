##
##  r k f 5 4 . R  Runge-Kutta-Fehlberg
##


rkf54 <- function(f, a, b, y0, tol = .Machine$double.eps^0.5,
                               control = list(), ...) {
    stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1)
    if (tol < .Machine$double.eps)  tol <- .Machine$double.eps

    # control list handling
    cntrl <- list(hmin = 0.001, hmax = 0.250, jmax = 200)
    nmsCo <- match.arg(names(control), choices = names(cntrl), several.ok = TRUE)
    if (!is.null(names(control))) cntrl[nmsCo] <- control

    fun <- match.fun(f)
    f   <- function(x, y) fun(x, y, ...)

    # coefficient matrices
    a2 <- 1/4; b2 <- 1/4
    a3 <- 3/8; b3 <- 3/32; c3 <- 9/32
    a4 <- 12/13; b4 <- 1932/2197; c4 <- -7200/2197; d4 <- 7296/2197
    a5 <- 1; b5 <- 439/216; c5 <- -8; d5 <- 3680/513; e5 <- -845/4104
    a6 <- 1/2; b6 <- -8/27; c6 <- 2; d6 <- -3544/2565; e6 <- 1859/4104; f6 <- -11/40

    r1 <- 1/360; r3 <- -128/4275; r4 <- -2197/75240; r5 <- 1/50; r6 <- 2/55
    n1 <- 25/216; n3 <- 1408/2565; n4 <- 2197/4104; n5 <- -1/5

    # Initialize solution vectors
    j <- 1
    T <- c(a)
    Y <- c(y0)

    # Initialize control parameters
    hmin <- cntrl$hmin
    hmax <- cntrl$hmax
    jmax <- cntrl$jmax
    h <- 0.8*hmin + 0.2*hmax
    br <- b - tol * abs(b)

    while (T[j] < b) {
        if (T[j] + h > br)  h <- b - T[j]
        tj <- T[j]; yj <- Y[j]

        k1 <- h * f(tj, yj)
        k2 <- h * f(tj + a2*h, yj + b2*k1)
        k3 <- h * f(tj + a3*h, yj + b3*k1 + c3*k2)
        k4 <- h * f(tj + a4*h, yj + b4*k1 + c4*k2 + d4*k3)
        k5 <- h * f(tj + a5*h, yj + b5*k1 + c5*k2 + d5*k3 + e5*k4)
        k6 <- h * f(tj + a6*h, yj + b6*k1 + c6*k2 + d6*k3 + e6*k4 + f6*k5)

        err  <-  abs(r1*k1 + r3*k3 + r4*k4 + r5*k5 + r6*k6)
        ynew <- yj + n1*k1 + n3*k3 + n4*k4 + n5*k5

        # Convergence condition
        if (err < tol || h < 2*hmin) {
            Y <- c(Y, ynew)
            if (tj + h > br) T <- c(T, b)
            else             T <- c(T, tj + h)
            j <- j+1;
        }

        # Compute next step length
        if (err == 0) s <- 0 else s <- s <- (tol*h/(2*err))^0.25
        if (s < 0.1) s <- 0.1
        if (s > 4.0) s <- 4.0

        h <- s * h
        if (h > hmax)  h <- hmax
        if (h < hmin)  h <- hmin

        if (j >= jmax) {
            warning(paste("Maximum number of steps reached:", j))
            break
        }
    }
    return(list(x = T, y = Y))
}
