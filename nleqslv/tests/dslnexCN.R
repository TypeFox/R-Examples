
library("nleqslv")

# Dennis & Schnabel,1996,"Numerical methods for unconstrained optimization and nonlinear equations", SIAM
# example 6.5.1 page 149

dslnex <- function(x) {
    y <- numeric(2)
    y[1] <- x[1]^2 + x[2]^2 - 2
    y[2] <- exp(x[1]-1) + x[2]^3 - 2
    y
}


xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart
fstart

do.print.xf <- TRUE

print.result <- function(z) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

sink("dslnexCN-num.txt")
for( z in c("dbldog","pwldog") ) {  # double dogleg, Powell (single) dogleg
    for( delta in c(-1.0, -2.0) ) { # Cauchy step , Newton step
        znlq <- nleqslv(xstart, dslnex, global=z, control=list(btol=.01,delta=delta, trace=1))
        print.result(znlq)
    }
}
sink()

sink("dslnexCN-char.txt")
for( z in c("dbldog","pwldog") ) {  # double dogleg, Powell (single) dogleg
    for( delta in c("cauchy", "newton") ) { # Cauchy step , Newton step
        znlq <- nleqslv(xstart, dslnex, global=z, control=list(btol=.01,delta=delta,trace=1))
        print.result(znlq)
    }
}
sink()

z1 <- readLines(con="dslnexCN-num.txt")
z2 <- readLines(con="dslnexCN-char.txt")

all.equal(z1,z2)
