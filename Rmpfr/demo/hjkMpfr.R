## some platforms hit zero exactly on the first step:
## if so the estimated precision is 2/3.
cyq.f <- function (x) {
    rv <- cyq.res(x)
    ##  mm <- length(rv)
    ##  f <- 0
    ##  for (ii in 1:mm) f <- f+rv[ii]*rv[ii]
    ## f <- sum(rv*rv)
    f <- crossprod(rv)
}

cyq200.f <- function (xx) {
    rv <- cyq200.res(xx)
    ##  mm <- length(rv)
    ##  f <- 0
    ##  for (ii in 1:mm) f <- f+rv[ii]*rv[ii]
    ## f <- sum(rv*rv)
    ##  f <- crossprod(rv)
    f <- sum(rv*rv)
}


cyq.res <- function (x) {
    ## Fletcher's chebyquad function m = n -- residuals
    n <- length(x)
    res <- rep(0,n)
    ##   res <- mpfrArray(rep(0,n), 200, dim=n) # initialize
    for (i in 1:n) { ## loop over resids
        rr <- 0
        for (k in 1:n) {
            z7 <- 1
            z2 <- 2*x[k]-1
            z8 <- z2
            j <- 1
            while (j<i) {
                z6 <- z7
                z7 <- z8
                z8 <- 2*z2*z7-z6 # recurrence to compute Chebyshev polynomial
                j <- j+1
            }                           # end recurrence loop
            rr <- rr+z8
        }                               # end loop on k
        rr <- rr/n
        if (2*trunc(i/2) == i) { rr <- rr + 1/(i*i - 1) }
        res[i] <- rr
    }                                   # end loop on i
    res
}

cyq200.res <- function (x) {
    ## Fletcher's chebyquad function m = n -- residuals
    n <- length(x)
    res <- mpfrArray(rep(0,n), 200, dim=n) # initialize
    for (i in 1:n) {                       #loop over resids
        rr <- 0
        for (k in 1:n) {
            z7 <- 1
            z2 <- 2*x[k]-1
            z8 <- z2
            j <- 1
            while (j<i) {
                z6 <- z7
                z7 <- z8
                z8 <- 2*z2*z7-z6 # recurrence to compute Chebyshev polynomial
                j <- j+1
            }                           # end recurrence loop
            rr <- rr+z8
        }                               # end loop on k
        rr <- rr/n
        if (2*trunc(i/2) == i) { rr <- rr + 1/(i*i - 1) }
        res[i] <- rr
    }                                   # end loop on i
    res
}

## JN need to make 1 an mpfr number of same length as x0 -- how to get precision
##	scale <- max(1, sqrt(sum(x0^2)))
##	scale <- max(mpfr(1, getPrec(x0)[1]), sqrt(sum(x0^2)))
cat("Minimize cyq for 2 parameters\n")
n <- 2
x <- 1:n
x <- x/(n+1) # Initial value suggested by Fletcher
##a2 <- hjkMpfr(x, cyq.f, control=list(trace=TRUE))
##print(a2)
x200 <- mpfrArray(1:n, 200)
x200 <- x200/mpfr(n+1,200)
cat("longstart:"); x200
cat("test residuals\n")
(rv0 <- cyq200.res(x200))

cat("test function\n")
(f0 <- cyq200.f(x200))

a2long <- hjkMpfr(x200, cyq200.f, control=list(info=TRUE))
print(a2long)

n <- 6
cat("Minimize cyq for ", n," parameters:\n")
x200 <- mpfr(1:n, 200)/(n+1)
cat("longstart:"); x200
cat("test residuals\n")
(rv0 <- cyq200.res(x200))

cat("test function\n")
(f0 <- cyq200.f(x200))

## Now this is *REALLY* very expensive {takes many minutes}:
a6long <- hjkMpfr(x200, cyq200.f, control=list(info=TRUE))
a6long
