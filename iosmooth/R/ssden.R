ssden <- function(x, y, bw) {
    dftRe <- function(s) {
		mean(y * cos(s * x))
	}
	dftIm <- function(s) {
		mean(y * sin(s * x))
    }
    n <- length(x)
    
    dftRe <- Vectorize(dftRe, "s")
    dftIm <- Vectorize(dftIm, "s")


	fhat <- function(xx) {
		fint <- function(s) kappaInDf(s, 1/bw)*(dftRe(s)*cos(s*xx) + dftIm(s)*sin(s*xx))
		integrate(fint, -2/bw, 2/bw)$value/(2*pi)
	}
	Vectorize(fhat, "xx")
}	
