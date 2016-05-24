# Broyden banded function

library("nleqslv")

brdban <- function(x) {
	ml <- 5
	mu <- 1
	n <- length(x)
    y <- numeric(n)

	for( k in 1:n ) {

		k1 <- max(1, k - ml)
		k2 <- min(n, k + mu)

		temp = 0.0
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

xsol <- c( -0.42830,  -0.47660,  -0.51965,  -0.55810,  -0.59251,
           -0.62450,  -0.62324,  -0.62139,  -0.62045,  -0.58647  )

fsol <- brdban(xsol)

znlq <- nleqslv(xstart, brdban, global="dbldog", method="Newton",
                control=list(trace=0,ftol=1e-8,xtol=1e-8,btol=1e-2,delta=-1.0))
znlq$termcd                 # should be 2 for x values within tolerance
all(abs(znlq$fvec)<=1e-7)   # may not have achieved ftol

xstart <- -2*rep(1,n)
znlq <- nleqslv(xstart, brdban, global="dbldog", method="Newton",
                control=list(trace=0,ftol=1e-8,xtol=1e-8,btol=1e-2,delta=-1.0))
znlq$termcd
all(abs(znlq$fvec)<=1e-8)

znlq <- nleqslv(xstart, brdban, global="dbldog",
                control=list(trace=0,ftol=1e-8,xtol=1e-8,btol=1e-2,delta=-1.0))
znlq$termcd                 # should be 2 for x values within tolerance
all(abs(znlq$fvec)<=1e-7)   # may not have achieved ftol

xstart <- -2*rep(1,n)
znlq <- nleqslv(xstart, brdban, global="dbldog",
                control=list(trace=0,ftol=1e-8,xtol=1e-8,btol=1e-2,delta=-1.0))
znlq$termcd
all(abs(znlq$fvec)<=1e-8)
