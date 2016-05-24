bwadap.ts <- function(x, Kn = 5, c.thresh = 2, ...) {
	n <- length(x)
	thresh <- c.thresh*sqrt(log(n, 10)/n)
	
	ac <- as.vector(acf(x, type="correlation", plot=FALSE, 
		lag.max = floor(n/2))$acf)
	l <- length(ac)
	
	pos <- 1
	
	while(pos < n/2)
	{
		npos <- match(TRUE, abs(ac[pos:l]) < thresh)
		if( is.na(npos) ) break;
		
		pos <- pos+npos-1
		if( pos+Kn-1 > floor(n/2) ) break;

		if(all(abs(ac[pos:(pos+Kn-1)]) < thresh)){
			return(pos-2)
		} else {
			npos <- match(FALSE, abs(ac[pos:(pos+Kn-1)]) < thresh)
			pos <- pos + npos - 1
		}
	}

	warning("No bandwidth found")		
	return(n/2)
}

bwadap.numeric <- function(x, smax=13.49/IQR(x), n.points = 1000, Kn = 1.349*5/IQR(x), c.thresh = 2, ...) {	
	n <- length(x)
	if(n <= 2) stop("x must have length greater than 2")

	thresh <- c.thresh*sqrt(log(n, 10)/n)
	
	dft <- function(s) {
		dftval <- complex( real = sum(cos(s * x)), imaginary = sum(sin(s * x)) )
		return(Mod(dftval)/n)
	}
	
	dft <- Vectorize(dft, "s")
	
	svals <- seq(0, smax, length.out=n.points)
	ftvals <- dft(svals)
	
	pos <- 1 
	
	while(pos < n.points)
	{
		npos <- match(TRUE, ftvals[pos:n.points] < thresh) + pos - 1
		if( is.na(npos) ) break;
		
		start.s <- uniroot(function(s) dft(s) - thresh, svals[c(npos-1,npos)])$root
		
		upcrosspos <- match(TRUE, ftvals[npos:n.points] > thresh)
		
		if( is.na(upcrosspos) ) {
			if(svals[n.points] - start.s < Kn) {
				warning("Not able to check Kn units beyond selected bandwidth")
			}
			return(1/start.s)
		}
		
		end.s <- uniroot(function(s) dft(s) - thresh, svals[c(npos + upcrosspos - 2, npos + upcrosspos - 1)])$root
		
		if( end.s - start.s > Kn ) {
			return(1/start.s)
		} else {
			pos <- npos + upcrosspos - 1
		}
	}

	warning("No bandwidth found")
	return(NA)
}

bwadap <- function(x, ...) UseMethod("bwadap")


