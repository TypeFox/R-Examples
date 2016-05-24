# Dennis & Schnabel,1996,"Numerical methods for unconstrained optimization and nonlinear equations", SIAM
# example 6.5.1 page 149

library(nleqslv)

dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}

jacdsln <- function(x) {
    n <- length(x)
    Df <- matrix(numeric(n*n),n,n)
    Df[1,1] <- 2*x[1]
    Df[1,2] <- 2*x[2]
    Df[2,1] <- exp(x[1]-1)
    Df[2,2] <- 3*x[2]^2

    Df
}

do.print.xf <- FALSE
do.trace <- 0

print.result <- function(z) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

xstart <- c(2,.5)

z <- nleqslv(xstart,dslnex, jacobian=TRUE, control=list(trace=do.trace))
print.result(z)
all.equal(z$jac,jacdsln(z$x), tolerance=0.05)

z <- nleqslv(xstart,dslnex,jacdsln, jacobian=TRUE, control=list(trace=do.trace))
print.result(z)
all.equal(z$jac,jacdsln(z$x), tolerance=0.05)

z <- nleqslv(xstart,dslnex, method="Newton", jacobian=TRUE, control=list(trace=do.trace))
print.result(z)
all.equal(z$jac,jacdsln(z$x), tolerance=10^3*.Machine$double.eps^0.5)

z <- nleqslv(xstart,dslnex, jacdsln, method="Newton", jacobian=TRUE, control=list(trace=do.trace))
print.result(z)
identical(z$jac,jacdsln(z$x))
