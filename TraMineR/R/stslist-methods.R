## ===========================
## Methods for stslist objects
## ===========================

print.stslist <- function(x,format='STS', extended=FALSE, ...) {
	if (format=='STS') {
		if (extended==FALSE) {
			void <- attr(x,"void")
			x <- seqconc(x, void=void)
			print(x, quote=FALSE, ...)
		} else NextMethod("print")
	}
	if (format=='SPS') {
		x <- seqconc(x, void=attr(x,"void"))
		
		if (extended==FALSE)
			x <- suppressMessages(seqformat(x,from='STS', to='SPS', compressed=TRUE))
		else if (extended==TRUE)
			x <- suppressMessages(seqformat(x,from='STS', to='SPS', compressed=FALSE))

		print(x, quote=FALSE)
	}
}

## plot.stslist <- function(x,...) {
##	seqiplot(x)
## }

"[.stslist" <- function(x,i,j,drop=FALSE) {
	## Specialized only for column subscript
	## If one column we keep the original data.frame method
	## Otherwise we copy attributes and update "start" value
	if (!missing(j) && length(j)>1) {
		## Storing the attributes
		x.attributes <- attributes(x)

		## Applying method
	     x <- NextMethod("[")
	
		## Adapting column names
		x.attributes$names <- x.attributes$names[j]

		## Redefining attributes
		attributes(x) <- x.attributes

	     attr(x,"start") <- x.attributes$start-1+j[1]

		if (!missing(i)) {
			attr(x,"row.names") <- attr(x,"row.names")[i]
			attr(x,"weights") <- attr(x,"weights")[i]
		}

		return(x)
	}

	x <- NextMethod("[")

	if (!missing(i))
		attr(x,"weights") <- attr(x,"weights")[i]
	
	return(x)
 }


## "[.stslist" <- function(x,...) {
## 	NextMethod("[")
## }

Math.stslist <- function(...){
 stop("Invalid operation on sequences")
}


