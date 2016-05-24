baseline.shirley <- function(spectra, t = NULL, limits = NULL, maxit = 50, err = 1e-6) {
# INPUT:
# spectra - matrix with y in first row (y = spectra[1,])
# t       - vector with x coordinates (OPTIONAL)
# limits  - spectrum limits
# maxit   - maximum number of iterations
# err     - cut-off error
#
# OUTPUT:
# baseline  - proposed baseline
# corrected - baseline corrected spectra

   	object <- spectra[1,]
	npt <- length(object)
    if ( is.null(t) ) t <- 1:npt
    
	if ( is.null(limits) ) {
		limits <- list( x = c(t[1], t[npt]), y = c(spectra[1,1], spectra[1,npt]) )
		}
	
	lowLim 	<- min(limits$y)
	BGND 	<- rep.int(lowLim, npt)
	SumBGND <- sum(BGND)
	SumRTF 	<- sum(object)
	RangeY 	<- diff(range(limits$y))
	
	if ( diff(limits$y) > 0 ) {
		nloop <- 0
		repeat {
			nloop <- nloop + 1
			for ( idx in rev(seq_len(npt))) {
				BGND[npt-idx+1] <- ( (RangeY/(SumRTF-sum(BGND)))*(sum(object[npt:idx])-sum(BGND[npt:idx])) ) + lowLim	
			}
			if ( ( abs( (sum(BGND)-SumBGND)/SumBGND ) < err ) || nloop > maxit ){ break }				
			SumBGND <- abs(sum(BGND))
		}

	}
	else {	
		nloop <- 0
		repeat {
			nloop <- nloop + 1
			for ( idx in seq_len(npt)) {
				BGND[idx] <- ( (RangeY/(SumRTF-sum(BGND)))*(sum(object[idx:npt])-sum(BGND[idx:npt])) ) + lowLim	
			}
			if ( ( abs( (sum(BGND)-SumBGND)/SumBGND ) < err ) || nloop > maxit ){ break }				
			SumBGND <- abs(sum(BGND))
		}
	}
		
	baseline <- matrix(data = BGND, nrow = 1)
	
	list(baseline = baseline, corrected = spectra - baseline)
}


