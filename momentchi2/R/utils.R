#A collection of methods used to compute moments/cumulants.

#ensure the coefficients are all positive
#returns TRUE if error, FALSE if no error
checkCoeffsArePositiveError <- function(coeff){
	if (length(coeff)==0){
        return(TRUE)
	}
	
	if ( !(all(coeff>0)) ){
        return(TRUE)
	}
    return(FALSE)
}

getCoeffError <- function(coeff){
	if (length(coeff)==0){
		return("empty coefficient vector.")
	}
	
	if ( !(all(coeff>0)) ){
		return("not all coefficients > 0.")
	}
    return("unknown error.")
}


#ensure the x-values are all positive
checkXvaluesArePositiveError <- function(x){
	if (length(x)==0){
        return(TRUE)
	}
	
	if ( !(all(x>0)) ){
        return(TRUE)
	}
    return(FALSE)
}

getXvaluesError <- function(x){
	if (length(x)==0){
		return("empty x vector.")
	}
	
	if ( !(all(x>0)) ){
		return("not all x-values > 0.")
	}
    return("unknown error.")
}

