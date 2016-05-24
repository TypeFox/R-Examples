## ===========================
## Methods for stsstatd objects
## ===========================

print.stslist.freq <- function(x, digits=2, width=1, ...) {
	table <- attr(x,"freq")
	print(table, digits=digits, width=width, ...)
}

"[.stslist.freq" <- function(x,i,j,drop=FALSE) {
	## Column subscript are not allowed
	if (!missing(j))
		stop(" [!] Column subscripts are not allowed", call.=FALSE)

	x <- NextMethod("[")

	if (!missing(i)) {
		attr(x,"weights") <- attr(x,"weights")[i]
		attr(x,"freq") <- attr(x,"freq")[i,]
	}
	
	return(x)
 }
