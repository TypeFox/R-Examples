##
##  e x p i n t . R  Exponential Integral
##


expint <- function(x) {
    stopifnot(is.numeric(x) || is.complex(x))
    eps <- .Machine$double.eps

    x <- c(x)
    n <- length(x)
    y <- numeric(n)

    p <- c(-3.602693626336023e-09, -4.819538452140960e-07, -2.569498322115933e-05,
           -6.973790859534190e-04, -1.019573529845792e-02, -7.811863559248197e-02,
           -3.012432892762715e-01, -7.773807325735529e-01,  8.267661952366478e+00)
    polyv <- polyval(p, Re(x))

    # series expansion
    k <- which(abs(Im(x)) <= polyv)
    if (length(k) != 0) {
       # initialization
       egamma <- 0.57721566490153286061
       xk <- x[k]
       yk <- -egamma - log(xk +0i)
       j <- 1
       pterm <- xk
       term  <- xk

       while (any(abs(term) > eps)) {
          yk <- yk + term
          j  <- j + 1
          pterm <- -xk * pterm / j
          term  <- pterm / j
       }
       y[k] <- yk
    }

    # continued fraction
    k <- which( abs(Im(x)) > polyv )
    if (length(k) != 0) {
        m <- 1  # we're calculating E1(x)

        # initialization
        xk <- x[k]
        nk <- length(xk)
        am2 <- numeric(nk)
        bm2 <- rep(1, nk)
        am1 <- rep(1, nk)
        bm1 <- xk;
        f <- am1 / bm1
        oldf <- rep(Inf, nk)
        j <- 2

        while (any(abs(f - oldf) > (100 * eps) * abs(f))) {
            alpha <- m - 1 + (j/2)

            # calculate the recursion coefficients
            a <- am1 + alpha * am2
            b <- bm1 + alpha * bm2

            # save new normalized variables for next pass
            am2 <- am1 / b
            bm2 <- bm1 / b
            am1 <- a / b
            bm1 <- 1
            f <- am1
            j <- j + 1

            # calculate the coefficients for j odd
            alpha <- (j-1)/2
            beta <- xk
            a <- beta * am1 + alpha * am2
            b <- beta * bm1 + alpha * bm2
            am2 <- am1 / b
            bm2 <- bm1 / b
            am1 <- a / b
            bm1 <- 1
            oldf <- f
            f <- am1
            j <- j+1
        }
        y[k] <- exp(-xk) * f - 1i * pi * ((Re(xk) < 0) & (Im(xk) == 0))
    }

    if (all(Im(y) == 0)) y <- Re(y)
    return(y)  
}

expint_E1 <- expint           # E1()


expint_Ei <- function(x) {    # Ei()
    stopifnot(is.numeric(x) || is.complex(x))
    # y <- -expint(-x) + sign(Im(x)) * pi * 1i
    y <- ifelse(sign(Im(x)) <= 0, -expint(-x) - pi*1i, -expint(-x) + pi*1i)
    if (all(Im(y) == 0)) y <- Re(y)
    return(y)
}

li <- function(x) {
    stopifnot(is.numeric(x) || is.complex(x))
    y <- expint_Ei(log(x + 0i))
    if (all(Im(y) == 0)) y <- Re(y)
    return(y)
}
