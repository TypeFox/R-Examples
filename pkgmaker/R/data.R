# Data related functions 
# 
# Author: Renaud Gaujoux
# Creation: 18 Jun 2012
###############################################################################

#' Generating Names
#' 
#' Generates names or dimnames for objects.
#' 
#' @param x object whose names are generated.
#' 
#' @export
#' @rdname addnames
addnames <- function(x, ...){
	UseMethod('addnames')
}

#' @S3method addnames default
#' @rdname addnames
addnames.default <- function(x, ...){
	if( is.vector(x) ) addnames.vector(x, ...)
	else 
		stop("no applicable method for 'addnames' applied to an object of class ", class(x))
}

#' @param prefix prefix string to use. A vector can be used to specify a prefix 
#' for each dimension of \code{x}. 
#' Names are build as \code{<prefix><sep><index>}.
#' @param sep separator used between the prefix and the numeric index. 
#' @param ... extra arguments to allow extension and passed to the next method.
#' 
#' @S3method addnames vector
#' @rdname addnames
addnames.vector <- function(x, prefix='x', sep='', ...){
	names(x) <- paste(prefix, 1:length(x), sep=sep) 
	x
} 


#' @S3method addnames array
#' @rdname addnames
addnames.array <- function(x, prefix=letters[1:length(dim(x))], sep='', ...){
	d <- dim(x)
	# recycle prefix if necessary
	if( length(prefix) != length(d) )
		prefix <- rep(prefix, length(d))
	
	dimnames(x) <- 
			lapply(seq_along(d), function(i) paste(prefix[i], 1:d[i], sep=sep))
	x
} 

#' @S3method addnames matrix
#' @rdname addnames
addnames.matrix <- function(x, prefix=c('row', 'col'), ...){
	addnames.array(x, prefix=prefix, ...)
}
