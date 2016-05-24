# equation.R -- version 2011-01-03
m <- 10000; n <- 10; trials <- 200 

X <- array(rnorm(m * n), dim = c(m, n))
y <- rnorm(m)

# QR
system.time({ 
    for (r in seq(trials)) sol1 <- qr.solve(X, y)
})

# form (X'X) and (X'y)
system.time({ 
    for (r in seq(trials)) 
        sol2 <- solve(crossprod(X), crossprod(X, y))
})

# cholesky
system.time({ 
    for (r in seq(trials)) {
        C <- chol(crossprod(X))
        rhs <- crossprod(X, y)
        sol3 <- backsolve(C, forwardsolve(t(C), rhs))
    }
})

# check
stopifnot(all.equal(as.numeric(sol1), as.numeric(sol2)))
stopifnot(all.equal(as.numeric(sol2), as.numeric(sol3)))