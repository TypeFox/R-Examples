
library("nleqslv")

# Trigonometric function
trig <- function(x) {
	n <- length(x)
	y <- cos(x)
	s <- sum(y)
	y <- n - s + c(1:n) * (1-y) - sin(x)

	y
}

trigjac <- function(x) {
	n <- length(x)
	J <- matrix(numeric(n*n),n,n)

	for (p in 1:n) {
		J[,p]  <- sin(x[p])
		J[p,p] <- (p+1) * sin(x[p]) - cos(x[p])
	}

	J
}

do.print.xf <- FALSE

print.result <- function(z) {
    if( do.print.xf ) {
        print(z$x)
        print(z$fvec)
    }
    print(z$message)
    print(all(abs(z$fvec)<=1e-8))
}

n <- 10
xstart <- rep(1,n)/n
fstart <- trig(xstart)

znlm <- nleqslv(xstart, trig, global="dbldog", control=list(trace=0))
print.result(znlm)
