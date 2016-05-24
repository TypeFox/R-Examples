
library(nleqslv)

# Dennis & Schnabel,1996,"Numerical methods for unconstrained optimization and nonlinear equations", SIAM
# example 6.5.1 page 149

dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
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

xcmp.result <- function(z1,z2)  all(abs(z1$x-z2$x) <= 1e-8)

xstart <- c(2,0.5)
hnlq1 <- nleqslv(xstart, dslnex, global="hook", control=list(btol=.01,delta="cauchy", trace=do.trace))
hnlq2 <- nleqslv(xstart, dslnex, global="hook", control=list(btol=.01,delta="newton", trace=do.trace))
print.result(hnlq1)
print.result(hnlq2)
xcmp.result(hnlq1,hnlq2)

dnlq1 <- nleqslv(xstart, dslnex, global="dbldog", control=list(btol=.01,delta="cauchy", trace=do.trace))
dnlq2 <- nleqslv(xstart, dslnex, global="dbldog", control=list(btol=.01,delta="newton", trace=do.trace))
print.result(dnlq1)
print.result(dnlq2)
xcmp.result(dnlq1,dnlq2)
xcmp.result(hnlq1,dnlq1)

xstart <- c(1.1,1.1)
hnlq1 <- nleqslv(xstart, dslnex, global="hook", control=list(btol=.01,delta="cauchy", trace=do.trace))
hnlq2 <- nleqslv(xstart, dslnex, global="hook", control=list(btol=.01,delta="newton", trace=do.trace))
print.result(hnlq1)
print.result(hnlq2)
xcmp.result(hnlq1,hnlq2)

dnlq1 <- nleqslv(xstart, dslnex, global="dbldog", control=list(btol=.01,delta="cauchy", trace=do.trace))
dnlq2 <- nleqslv(xstart, dslnex, global="dbldog", control=list(btol=.01,delta="newton", trace=do.trace))
print.result(dnlq1)
print.result(dnlq2)
xcmp.result(dnlq1,dnlq2)
xcmp.result(hnlq1,dnlq1)
