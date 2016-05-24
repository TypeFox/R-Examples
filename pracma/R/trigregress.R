##
##  t r i g r e g r e s s . R  Trigonometric Regression
##


trigPoly <- function(x, m) {
    stopifnot(is.numeric(x), is.numeric(m), length(m) == 1, m >= 0)
    if (m == 0)
        return(list(a0 = mean(x), a = c(), b = c()))

    n <- length(x)
    t <- seq(0, 2*(n-1)*pi/n, length.out = n)
    t <- as.matrix(t)

    a <- b <- numeric(m)
    for (j in 1:m) {
        a[j] <- x %*% cos(j*t)
        b[j] <- x %*% sin(j*t)
    }
    a <- 2*a/n
    b <- 2*b/n
    a0 <- sum(x)/n
    if (n == 2*m) a[m] <- a[m]/2

    return(list(a0 = a0, a = a, b = b))
}


trigApprox <- function(t, x, m) {
    stopifnot(is.numeric(t))

    tP <- trigPoly(x, m)
    a0 <- tP$a0
    a <- tP$a; b <- tP$b

    y <- a0
    for (j in 1:m) {
        y <- y + a[j]*cos(j*t) + b[j]*sin(j*t)
    }
    return(y)
}
