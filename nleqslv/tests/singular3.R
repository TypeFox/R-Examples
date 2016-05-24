
library(nleqslv)

print.result <- function(z, do.print.xf=FALSE) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

# Powell cautionary example
# M.J.D. Powell, "A Hybrid Method for Nonlinear Equations",
# in Numerical methods for Nonlinear Algebraic Equations, ed. P. Rabinowitz, 1970.


f <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]
    y[2] <- 10*x[1]/(x[1]+.1) + 2*x[2]^2

    y
}

jac <- function(x) {
    fjac <- matrix(0,nrow=2,ncol=2)

    fjac[1, 1] <- 1
    fjac[1, 2] <- 0
    fjac[2, 1] <- 1/(x[1]+0.1)^2
    fjac[2, 2] <- 4*x[2]

    fjac
}

xstart <- c(3,1)
z1 <- nleqslv(xstart,f, method="Newton",control=list(trace=0,allowSingular=TRUE))
print.result(z1)
xstart <- c(3,0) # singular start
z2 <- nleqslv(xstart,f, method="Newton",control=list(trace=0,allowSingular=TRUE))
print.result(z2)
