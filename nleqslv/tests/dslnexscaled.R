# Dennis Schnabel example

library("nleqslv")

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

xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart
fstart

# a solution is c(1,1)

do.print.xf <- FALSE

print.result <- function(z) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

# Use our fixed scaling

znlq <- nleqslv(xstart, dslnex, jacdsln, global="dbldog", control=list(btol=.01,delta=-1.0,chkjac=TRUE,scalex=c(2,3)))
if(znlq$termcd == -10) stop("Internal error in check analytical jacobian")

# Broyden analytical jacobian
for( z in c("dbldog","pwldog") ) {  # double dogleg, Powell (single) dogleg
    for( delta in c(-1.0, -2.0) ) { # Cauchy step , Newton step
        znlq <- nleqslv(xstart, dslnex, jacdsln, global=z, control=list(btol=.01,delta=delta,chkjac=TRUE,scalex=c(2,3)))
        print.result(znlq)
    }
}

# Newton analytical jacobian
for( z in c("dbldog","pwldog") ) {  # double dogleg, Powell (single) dogleg
    for( delta in c(-1.0, -2.0) ) { # Cauchy step , Newton step
        znlq <- nleqslv(xstart, dslnex, jacdsln, method="Newton", global=z, control=list(btol=.01,delta=delta,chkjac=TRUE,scalex=c(2,3)))
        print.result(znlq)
    }
}
