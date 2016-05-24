
library(nleqslv)

# Brown almost linear function

brown <- function(x) {
  	n <- length(x)
	y <- numeric(n)

  	y[1:(n-1)] <- x[1:(n-1)] + sum(x[1:n]) - (n + 1)
  	y[n] <- prod(x[1:n]) - 1.0

	y
}

brownjac <- function(x) {
    n <- length(x)
    J <- matrix(1,nrow=n,ncol=n)
    diag(J) <- 2
	xprod <- prod(x)
    J[n,] <- xprod/x  # exact
    J
}

print.result <- function(z, do.print.xf=FALSE) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

for( m in c("Newton","Broyden") ) {
    for( n in c(50,100) ) {
        xstart <- rep(1,n)/2
        z <- nleqslv(xstart, brown, brownjac, method="Newton",
                        control=list(trace=0,ftol=1e-10,delta="cauchy",allowSingular=TRUE))
        print.result(z)
    }
}
