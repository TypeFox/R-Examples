

setMethod("show", "cprobd", function(object) {

	cat(" [>] context:", object@context, "\n")
	show(object@.Data)
}
)


setMethod("round", "cprobd", function(x, digits=0) {

	x@.Data <- round(x@.Data, digits=digits)
	return(x)
}
)

setMethod("show", "cprobd.list", function(object) {

	for (i in 1:length(object)) {
		show(object[[i]])
		if (i < length(object)) { cat("\n") }
	}
}
)


setMethod("[", "cprobd.list", function(x, i, drop = TRUE) {
	
		A <- x@alphabet
		cpal <- x@cpal
		labels <- x@labels
		x <- x@.Data
		
		res <- new("cprobd.list", x[i], alphabet=A, cpal=cpal, labels=labels)
		return(res)
}
)

setMethod("plot", "cprobd.list", function(x, ...) {
	seqdata <- unlist(lapply(x, function(x) { x@context }))
	seqdata <- seqdef(seqdata, alphabet=x@alphabet, cpal=x@cpal, stlab=x@labels, nr="#")
	plot(seqdata, ...)
}
)
 
