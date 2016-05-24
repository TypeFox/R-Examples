## ============================
## Methods for stsmeant objects
## ============================

print.stslist.meant <- function(x, digits=2, ...) {

	cn <- colnames(x)
	## Conversion for printing without attributes
	x <- as.matrix(x[1:nrow(x),1:ncol(x)])
	colnames(x) <- cn

	NextMethod("print", digits=digits,...)
}
