sv <- function(the.vect, rownum=1){ #select variable
	#avoid going beyond the vector.
	#If we have only one entry, we may
	# assume that this is valid for all
	if(length(the.vect)==1) rownum=1

	stopifnot(is.null(the.vect) || rownum <= length(the.vect))

	if(is.null(the.vect)){
		to.ret <- "NULL" 
	} else if(is.na(the.vect[rownum])){
		to.ret <- "NULL" 
	} else {
		to.ret <- the.vect[rownum]
	}
	return(to.ret)
}

