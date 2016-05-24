## ===========================
## Methods for stsstatd objects
## ===========================

print.stslist.modst <- function(x, digits=2, ...) {

	cat("  [Modal state sequence]\n")
	## x <- NextMethod("print",...)
	print.stslist(x,...)	

	cat("\n  [State frequencies]\n")
	print(attr(x,"Frequencies"), digits=digits)
}

"[.stslist.modst" <- function(x,i,j,drop=FALSE) {
	## Specialized only for column subscript
	## If one column we keep the original data.frame method
	## Otherwise we copy attributes and update "start" value
	if (!missing(i))
		stop("row subscripts not allowed!", call.=FALSE)
	
	if (!missing(j)) {
		freq <- attr(x,"Frequencies")

		## Applying method
	     x <- NextMethod("[")
	
		## Adapting frequencies
		attr(x,"Frequencies") <- freq[j]
	}
	
	return(x)
 }
