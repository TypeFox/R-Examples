# Broyden banded function

library(nleqslv)

brdtri <- function(x) {
	n <- length(x)
    y <- numeric(n)

    # y[1] <- (3-2*x[1])*x[1] - 2*x[2] + 1
    # y[n] <- (3-2*x[n])*x[n] - x[n-1] + 1
    #
    # k <- 2:(n-1)
    # y[k] <- (3-2*x[k])*x[k] - x[k-1] - 2 * x[k+1] + 1
    #
    y <- (3-2*x)*x - c(0,x[-n]) - 2*c(x[-1],0) + 1
	y
}

brdtrijac <- function(x) {
    n <- length(x)
    J <- diag(3-4*x,n,n)
    J[row(J)==col(J)+1] <- -1
    J[row(J)==col(J)-1] <- -2
    J
}

options(echo=TRUE)

n <- 10
xstart <- -rep(1,n)
fstart <- brdtri(xstart)

z0 <- nleqslv(xstart,brdtri, method="Newton", global="dbldog")
z0$message

z1 <- nleqslv(xstart,brdtri, brdtrijac, method="Newton", global="dbldog",control=list(trace=0))
z1$message
all.equal(z1$x,z0$x)

z2 <- nleqslv(xstart,brdtri, brdtrijac, method="Newton", global="dbldog",control=list(trace=0,chkjac=TRUE))
z2$message
all.equal(z2$x,z0$x)

z3 <- nleqslv(xstart,brdtri, brdtrijac, method="Newton", global="dbldog",control=list(trace=0,dsub=1,dsuper=1))
z3$message
all.equal(z2$x,z0$x)

z4 <- nleqslv(xstart,brdtri, brdtrijac, method="Newton", global="dbldog",control=list(trace=0,dsub=1,dsuper=1,chkjac=TRUE))
z4$message
all.equal(z2$x,z0$x)
