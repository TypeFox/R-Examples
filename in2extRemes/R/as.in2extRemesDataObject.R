as.in2extRemesDataObject <- function(x) {

    out <- list()

    if(is.matrix( x) || is.data.frame(x)) {

	if(!is.data.frame(x)) x <- as.data.frame(x)
	out$data <- x

	if(is.null(colnames(x))) colnames(out$data) <- paste("x", 1:length(x), sep="")
	else if(is.null(colnames(out$data))) colnames(out$data) <- colnames(x)

    } else if(is.vector(x)) {

	out$data <- cbind(1:length(x), x)

	colnames(out$data) <- c("obs", deparse(substitute(x)))

    } else stop("as.in2extRemesDataObject: x must be a data frame, matrix or vector.")

    class(out) <- "in2extRemesDataObject"

    return(out)

} # end of 'as.in2extRemesDataObject' function.
