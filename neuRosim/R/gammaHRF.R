gammaHRF <-
function(x, FWHM=4, verbose=TRUE){

	th <- 0.242*FWHM
	
	1/(th*factorial(3)) * (x/th)^3 * exp(-x/th)
}

