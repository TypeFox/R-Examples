##
##  f z e r o . R
##


fzero <- function(f, x, ..., maxiter = 100, tol = .Machine$double.eps^(1/2)) {
    if (!is.numeric(x) || length(x) > 2)
        stop("Argument 'x' must be a scalar or a vector of length 2.")

    err <- try(fun <- match.fun(f), silent = TRUE)
    if (class(err) == "try-error") {
        stop("Argument function 'f' not known in parent environment.")
    } else {
        f <- function(x) fun(x, ...)
    }

    zin <- NULL
    if (length(x) == 2) {
        if (x[1] <= x[2]) {
            a <- x[1]; b <- x[2]
        } else {
            warning("Left endpoint bigger than right one: exchanged points.")
            a <- x[2]; b <- x[1]
        }
        zin <- .zeroin(f, a, b, maxiter = maxiter, tol = tol)

    } else {  # try to get b
        a <- x; fa <- f(a)
        if (fa == 0) return(list(x = a, fval = fa))
        if (a == 0) {
            aa <- 1
        } else {
            aa <- a
        }
        bb <- c(0.9*aa, 1.1*aa, aa-1, aa+1, 0.5*aa, 1.5*aa,
                -aa, 2*aa, -10*aa, 10*aa)
        for (b in bb) {
            fb <- f(b)
            if (fb == 0) return(list(x = b, fval = fb))
            if (sign(fa) * sign(fb) < 0) {
                zin <- .zeroin(f, a, b, maxiter = maxiter, tol = tol)
                break
            }
        }
    }

    if (is.null(zin)) {
        warning("No interval w/ function 'f' changing sign was found.")
        return(list(x = NA, fval = NA))
    } else {
        x1 <- zin$bra[1]; x2 <- zin$bra[2]
        f1 <- zin$ket[1]; f2 <- zin$ket[2]
        x0 <- sum(zin$bra)/2; f0 <- f(x0)
        if (f0 < f1 && f0 < f2) {
            return(list(x = x0, fval = f0))
        } else if (f1 <= f2) {
            return(list(x = x1, fval = f1))
        } else {
            return(list(x = x2, fval = f2))
        }
    }
}


.zeroin <- function(f, a, b, maxiter = 100, tol = 1e-07) {
    stopifnot(is.numeric(a), is.numeric(b),
                length(a) == 1, length(b) == 1)

    mu <- 0.5
    eps <- .Machine$double.eps
    info <- 0

    ## Set up and prepare ...
    if (a > b) {
        tt <- a; a <- b; b <- tt
    }
    fa <- f(a); fb <- f(b) 
    nfev <- 2
    if (sign(fa) * sign(fb) > 0)
        stop("Function must differ in sign at interval endpoints.")

    itype <- 1
    if (abs(fa) < abs(fb)) {
        u <- a; fu <- fa
    } else {
        u <- b; fu <- fb
    }
    d <- e <- u
    fd <- fe <- fu

    mba <- mu * (b - a)
    niter <- 0
    while (niter < maxiter) {
        if (itype == 1) {
            # The initial test
            if ((b-a) <= 2*(2*abs(u)*eps + tol)) {
                x <- u; fval <- fu
                info <- 1
                break
            }
            if (abs (fa) <= 1e3*abs(fb) && abs(fb) <= 1e3*abs(fa)) {
                # Secant step.
                c <- u - (a - b) / (fa - fb) * fu
            } else {
                # Bisection step.
                c <- 0.5 * (a + b)
            }
            d <- u; fd <- fu
            itype <- 5

        } else if (itype == 2 || itype == 3) {
            l <- length(unique(c(fa, fb, fd, fe)))
            if (l == 4) {
                # Inverse cubic interpolation.
                q11 <- (d - e) * fd / (fe - fd)
                q21 <- (b - d) * fb / (fd - fb)
                q31 <- (a - b) * fa / (fb - fa)
                d21 <- (b - d) * fd / (fd - fb)
                d31 <- (a - b) * fb / (fb - fa)
                q22 <- (d21 - q11) * fb / (fe - fb)
                q32 <- (d31 - q21) * fa / (fd - fa)
                d32 <- (d31 - q21) * fd / (fd - fa)
                q33 <- (d32 - q22) * fa / (fe - fa)
                c <- a + q31 + q32 + q33;
            }
            if (l < 4 || sign(c - a) * sign(c - b) > 0) {
                # Quadratic interpolation + newton
                a0 <- fa
                a1 <- (fb - fa)/(b - a)
                a2 <- ((fd - fb)/(d - b) - a1) / (d - a)
                # Modification 1: this is simpler and does not seem to be worse.
                c <- a - a0/a1
                if (a2 != 0) {
                    c <- a - a0/a1
                    for (i in 1:itype) {
                        pc <- a0 + (a1 + a2*(c - b))*(c - a)
                        pdc <- a1 + a2*(2*c - a - b)
                        if (pdc == 0) {
                            c <- a - a0/a1
                            break
                        }
                        c <- c - pc/pdc
                    }
                }
            }
            itype <- itype + 1 
            
        } else if (itype == 4) {
            # Double secant step.
            c <- u - 2*(b - a)/(fb - fa)*fu
            # Bisect if too far.
            if (abs (c - u) > 0.5*(b - a)) {
                c <- 0.5 * (b + a)
            }
            itype <- 5

        } else if (itype == 5) {
            # Bisection step.
            c <- 0.5 * (b + a)
            itype <- 2
        }

        # Don't let c come too close to a or b.
        delta <- 2 * 0.7 * (2 * abs(u) * eps + tol)
        if ((b - a) <= 2*delta) {
            c <- (a + b)/2
        } else {
            c <- max(a + delta, min(b - delta, c))
        }

        # Calculate new point.
        x <- c;
        fval <- fc <- f(c)
        niter <- niter + 1; nfev <- nfev + 1

        # Mod2: skip inverse cubic interpolation if nonmonotonicity is detected.
        if (sign(fc - fa) * sign(fc - fb) >= 0) {
          ## The new point broke monotonicity. 
          ## Disable inverse cubic.
          fe <- fc
        } else {
          e <- d; fe <- fd
        }

        # Bracketing.
        if (sign(fa) * sign(fc) < 0) {
            d <- b; fd <- fb
            b <- c; fb <- fc
        } else if (sign(fb) * sign(fc) < 0) {
            d <- a; fd <- fa
            a <- c; fa <- fc
        } else if (fc == 0) {
            a <- b <- c; fa <- fb <- fc
            info <- 1
            break
        } else {
            # This should never happen.
            stop("zeroin: zero point could not be bracketed")
        }

        if (abs(fa) < abs(fb)) {
          u <- a; fu <- fa
        } else {
          u <- b; fu <- fb
        }
        if (b - a <= 2*(2 * abs (u) * eps + tol)) {
          info <- 1
          break
        }

        # Skip bisection step if successful reduction.
        if (itype == 5 && (b - a) <= mba) {
          itype <- 2
        }
        if (itype == 2) {
          mba <- mu * (b - a)
        }
    } # endwhile

    return(list(bra = c(a, b), ket = c(fa, fb), niter = niter, info = info))
}
