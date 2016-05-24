# Function from package 'fda' (c) 2014

is.eqbasis <- function(basisobj1, basisobj2) {
	
	#  tests to see of two basis objects are identical
	
	eqwrd <- TRUE
	
	#  test type
	
	if (basisobj1$type != basisobj2$type) {
		eqwrd <- FALSE
		return(eqwrd)
	}
	
	#  test range
	
	if (any(basisobj1$rangeval != basisobj2$rangeval)) {
		eqwrd <- FALSE
		return(eqwrd)
	}
	
	#  test nbasis
	
	if (basisobj1$nbasis != basisobj2$nbasis) {
		eqwrd <- FALSE
		return(eqwrd)
	}
	
	#  test params
	
	if (any(basisobj1$params != basisobj2$params)) {
		eqwrd <- FALSE
		return(eqwrd)
	}

   #  test dropind
	
	if (any(basisobj1$dropind != basisobj2$dropind)) {
		eqwrd <- FALSE
		return(eqwrd)
	}
	
	return(eqwrd)
	
}
