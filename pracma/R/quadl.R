##
##  q u a d l . R  Adaptive Simpson Quadrature
##


quadl <- function(f, xa, xb, tol = .Machine$double.eps^0.5, trace = FALSE, ...)
{
    stopifnot(is.numeric(xa), length(xa) == 1, is.finite(xa),
              is.numeric(xb), length(xb) == 1, is.finite(xb))

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    if (xa == xb) return(xb-xa)
    else if (xa > xb) {
        tmp <- xa; xa <- xb; xb <- tmp
        rev_p <- TRUE
    } else
        rev_p <- FALSE

    eps <- .Machine$double.eps
    if (!is.finite(f(xa))) xa <- xa + 2*eps
    if (!is.finite(f(xb))) xb <- xb - 2*eps

    Q <- .adaptlob(f, xa, xb, tol, trace)
    if (rev_p) Q <- -1 * Q
    return(Q)
}


.adaptlob <- function(f, a, b, tol = tol, trace = trace)
{
    m <- (a+b)/2
    h <- (b-a)/2
    alpha <- sqrt(2/3)
    beta <- 1/sqrt(5)

    x1 <- 0.942882415695480
    x2 <- 0.641853342345781
    x3 <- 0.236383199662150
    x <- c(a, m-x1*h, m-alpha*h, m-x2*h, m-beta*h, m-x3*h,  m,
              m+x3*h, m+beta*h,  m+x2*h, m+alpha*h, m+x1*h, b)

    y <- f(x)
    fa <- y[1]
    fb <- y[13]
    i2 <- (h/6) * (y[1] + y[13] + 5*(y[5]+y[9]))
    i1 <- (h/1470) * (77*(y[1]+y[13]) + 432*(y[3]+y[11]) + 
                      625*(y[5]+y[9]) + 672*y[7])
    ab  <- h * (0.0158271919734802 * (y[1]+y[13]) +
                0.0942738402188500 * (y[2]+y[12]) +
                0.155071987336585  * (y[3]+y[11]) +
                0.188821573960182  * (y[4]+y[10]) +
                0.199773405226859  * (y[5]+y[9])  +
                0.224926465333340  * (y[6]+y[8])  +
                0.242611071901408  *  y[7])

    s <- sign(ab)
    if (s == 0) s <- 1
    erri1 <- abs(i1-ab)
    erri2 <- abs(i2-ab)
    R <- 1
    if (erri2 != 0) R <- erri1/erri2
    if (R > 0 && R < 1) tol <- tol/R
    ab <- s * abs(ab) * tol/.Machine$double.eps
    if (ab == 0) ab <- b-a

    Q <- .adaptlobstp(f, a, b, fa, fb, ab, trace)
}


.adaptlobstp <- function(f, a, b, fa, fb, ab, trace)
{
    h <- (b-a)/2
    m <- (a+b)/2
    alpha <- sqrt(2/3)
    beta <- 1/sqrt(5)
    mll <- m - alpha*h
    ml  <- m - beta*h
    mr  <- m + beta*h
    mrr <- m + alpha*h

    x <- c(mll, ml, m, mr, mrr)
    y <- f(x)

    fmll <- y[1]
    fml  <- y[2]
    fm   <- y[3]
    fmr  <- y[4]
    fmrr <- y[5]
    i2 <- (h/6) * (fa + fb + 5*(fml+fmr))
    i1 <- (h/1470) * (77*(fa+fb) + 432*(fmll+fmrr) + 625*(fml+fmr) + 672*fm)

    if ( ab+(i1-i2) == ab | mll <= a | b <= mrr ) {
        if ( (m <= a || b <= m) )
            warning("Required tolerance may not be met.")

        Q <- i1
        if (trace) cat(a, b-a, Q, "\n")

    } else {
        Q <- .adaptlobstp(f, a, mll, fa, fmll, ab, trace)   +
             .adaptlobstp(f, mll, ml, fmll, fml, ab, trace) +
             .adaptlobstp(f, ml, m, fml, fm, ab, trace)     +
             .adaptlobstp(f, m, mr, fm, fmr, ab, trace)     +
             .adaptlobstp(f, mr, mrr, fmr, fmrr, ab, trace) +
             .adaptlobstp(f, mrr, b, fmrr, fb, ab, trace)
    }
    return(Q)
}
