mscale <- function(u, na.rm=FALSE)
#	Scale M-estimator with 50% breakdown
#	Yohai (1987) Annals, Stromberg (1993) JASA.
#
#	GKS  2 June 1999
#	Revised 17 April 2010
{
	isna <- is.na(u)
	if(any(isna)) {
		if(na.rm) {
			if(any(!isna))
				u <- u[!isna]
			else
				return(NA)
		} else {
			return(NA)
		}	
	}
	if(mean(u==0) >= 0.5) return(0)
	U <- abs(u)
	s <- median(U)/0.6744898
	iter <- 0
	repeat {
		iter <- iter+1
		z <- u/0.212/s
		d1 <- mean(.rho.hampel(z))-3.75
		d2 <- mean(z*.psi.hampel(z))
		s <- s*(1+d1/d2)
		if(iter > 50) {
			warning("Max iterations exceeded")
			break
		}
		if(abs(d1/d2) < 1e-13) break
	}
	s	
}

.rho.hampel <- function(u, a = 1.5, b = 3.5, c = 8)
{
#	Integral of Hampel's redescending psi function (Hampel, Ronchetti,
#	Rousseeuw and Stahel, 1986, Robust Statistics, Wiley, page 150).
#	Default values are as in Stromberg (1993) JASA.
#
#	GKS  31 May 99
#
	U <- abs(u)
	A <- (U <= a)	#increasing
	B <- (U > a) & (U <= b)	#flat
	C <- (U > b) & (U <= c)	#descending
	D <- (U > c)	# zero
	rho <- U
	rho[A] <- (U[A] * U[A])/2
	rho[B] <- a * (U[B] - a/2)
	rho[C] <- a * (b - a/2) + a * (U[C] - b) * (1 - (U[C] - b)/(c - b)/2)
	rho[D] <- (a * (b - a + c))/2
	rho
}

.psi.hampel <- function(u, a = 1.5, b = 3.5, c = 8)
{
#	Hampel's redescending psi function (Hampel, Ronchetti,
#	Rousseeuw and Stahel, 1986, Robust Statistics, Wiley, page 150).
#	Default values are as in Stromberg (1993) JASA.
#
#	GKS  2 June 99
#
	U <- abs(u)
	B <- (U > a) & (U <= b)	#flat
	C <- (U > b) & (U <= c)	#descending
	D <- (U > c)	# zero
	psi <- u
	psi[B] <- sign(u[B]) * a
	psi[C] <- sign(u[C]) * a * (c - U[C])/(c - b)
	psi[D] <- 0
	psi
}