mytest_iterator <- function(testdata, FUN=c("sum", "prod", "sumpower"), OUT, MAR=c(1, 2), POW=0)
{
	
	allowed = c("sum", "prod", "sumpower")
	if (length(which(allowed == FUN)) == 0)
		stop("Function should be one of sum, prod or sumpower")
	
	allowed = c(1, 2)
	if (length(which(allowed == MAR)) == 0)
		stop("Margin should be 1 (row-wise) or 2 (column-wise)")
	
	if (is(testdata,"databel")) {
	
		if (FUN == "sum" || FUN == "prod") {
			out <- .Call("iterator", testdata@data, as.integer(0), as.integer(0), as.character(FUN), 
						as.character(OUT), as.integer(MAR), as.integer(0), package="DatABEL")
			return(out)
		}
		
		if (FUN == "sumpower") {
			if (!is.numeric(POW))
				stop("Exponent should be a numeric value")
			out <- .Call("iterator", testdata@data, as.integer(0), as.integer(0), as.character(FUN),
						as.character(OUT), as.integer(MAR), as.integer(1), as.double(POW), 
						package="DatABEL")
			return(out)
		}
	
	} else if (is(testdata,"gwaa.data")) {
		
		nids <- testdata@gtdata@nids
		nsnps <- testdata@gtdata@nsnps
		
		if (FUN == "sum" || FUN == "prod") {
			out <- .Call("iterator", as.raw(testdata@gtdata@gtps), as.integer(nids), as.integer(nsnps),
						as.character(FUN), as.character(OUT), as.integer(MAR), as.integer(0), 
						package="DatABEL")
			return(out)
		}
		
		if (FUN == "sumpower") {
			if (!is.numeric(POW))
				stop("Exponent should be a numeric value")
			out <- .Call("iterator", as.raw(testdata@gtdata@gtps), as.integer(nids), as.integer(nsnps),
						as.character(FUN), as.character(OUT), as.integer(MAR), as.integer(1),
						as.double(POW), package="DatABEL")
			return(out)
		}
		
	} else {
		stop("Data argument should be of DatABEL or old gwaa.data class")
	}
	
}