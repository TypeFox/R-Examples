# Broyden banded

library("nleqslv")

brdban <- function(x,ml=5,mu=1) {
    n <- length(x)
    y <- numeric(n)

    for( k in 1:n ) {

        k1 <- max(1, k - ml)
        k2 <- min(n, k + mu)

        temp <- 0.0
        for(j in k1:k2) {
            if ( j != k ) {
                temp <- temp + x[j] * (1.0 + x[j])
            }
        }

        y[k] <- x[k] * (2.0 + 5.0 * x[k]**2) + 1.0 - temp

    }
    y
}

n <- 10
xstart <- -rep(1,n)
ztol <- 1000*.Machine$double.eps

z1 <- nleqslv(xstart,brdban, method="Newton")
z2 <- nleqslv(xstart,brdban, method="Newton", control=list(dsub=5,dsuper=1))

cat("z1 termcd=",z1$termcd, "jcnt,fcnt=",z1$njcnt,z1$nfcnt,"\n")
cat("z2 termcd=",z2$termcd, "jcnt,fcnt=",z2$njcnt,z2$nfcnt,"\n")
z1$message
z2$message
all.equal(z2$x,z1$x)
all.equal(z2$x,z1$x, tolerance=ztol)

z1 <- nleqslv(xstart,brdban, ml=2,mu=2, method="Newton")
z2 <- nleqslv(xstart,brdban, ml=2,mu=2, method="Newton", control=list(dsub=2,dsuper=2))

cat("z1 termcd=",z1$termcd, "jcnt,fcnt=",z1$njcnt,z1$nfcnt,"\n")
cat("z2 termcd=",z2$termcd, "jcnt,fcnt=",z2$njcnt,z2$nfcnt,"\n")
z1$message
z2$message
all.equal(z2$x,z1$x, tolerance=ztol)

z3 <- nleqslv(xstart,brdban, ml=2,mu=2, method="Broyden")
z4 <- nleqslv(xstart,brdban, ml=2,mu=2, method="Broyden", control=list(dsub=2,dsuper=2))

cat("z3 termcd=",z1$termcd, "jcnt,fcnt=",z3$njcnt,z3$nfcnt,"\n")
cat("z4 termcd=",z2$termcd, "jcnt,fcnt=",z4$njcnt,z4$nfcnt,"\n")
z3$message
z4$message
all.equal(z3$x,z1$x)
all.equal(z4$x,z1$x)
all.equal(z4$x,z3$x, tolerance=ztol)
