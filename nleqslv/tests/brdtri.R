# Broyden tridiagonal

library("nleqslv")

brdtri <- function(x) {
    n <- length(x)
    y <- numeric(n)

    y[1] <- (3-2*x[1])*x[1] - 2*x[2] + 1
    y[n] <- (3-2*x[n])*x[n] - x[n-1] + 1

    k <- 2:(n-1)
    y[k] <- (3-2*x[k])*x[k] - x[k-1] - 2 * x[k+1] + 1

    y
}

n <- 100
xstart <- -rep(1,n)
ztol <- 1000*.Machine$double.eps

z1 <- nleqslv(xstart,brdtri, method="Newton")
z2 <- nleqslv(xstart,brdtri, method="Newton", control=list(dsub=1,dsuper=1))

cat("z1 termcd=",z1$termcd, "jcnt,fcnt=",z1$njcnt,z1$nfcnt,"\n")
cat("z2 termcd=",z2$termcd, "jcnt,fcnt=",z2$njcnt,z2$nfcnt,"\n")
z1$message
z2$message
all.equal(z2$x,z1$x)
all.equal(z2$x,z1$x, tolerance=ztol)

z1 <- nleqslv(xstart,brdtri, method="Newton")
z2 <- nleqslv(xstart,brdtri, method="Newton", control=list(dsub=1,dsuper=1))

cat("z1 termcd=",z1$termcd, "jcnt,fcnt=",z1$njcnt,z1$nfcnt,"\n")
cat("z2 termcd=",z2$termcd, "jcnt,fcnt=",z2$njcnt,z2$nfcnt,"\n")
z1$message
z2$message
all.equal(z2$x,z1$x, tolerance=ztol)

z3 <- nleqslv(xstart,brdtri, method="Broyden")
z4 <- nleqslv(xstart,brdtri, method="Broyden", control=list(dsub=1,dsuper=1))

cat("z3 termcd=",z1$termcd, "jcnt,fcnt=",z3$njcnt,z3$nfcnt,"\n")
cat("z4 termcd=",z2$termcd, "jcnt,fcnt=",z4$njcnt,z4$nfcnt,"\n")
z3$message
z4$message
all.equal(z3$x,z1$x)
all.equal(z4$x,z1$x)
all.equal(z4$x,z3$x, tolerance=ztol)
