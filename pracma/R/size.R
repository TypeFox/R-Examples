##
##  s i z e . R  Matlab size, numel, ndims, and isempty functions
##


size <- function(x, k) {
	if (length(x) == 0)    sz <- 0
	else if (is.vector(x)) sz <- c(1, length(x))
	else if (is.array(x))  sz <- dim(x)
	else                   sz <- NULL

	if (!missing(k)) {
		if (k > length(sz)) sz <- 1
		else if (k >= 1)    sz <- sz[k]
		else
			stop("Requested dimension 'k' is out of range.")
	}
	return(sz)
}

numel <- function(x) {
	sz <- size(x)
	if (!is.null(sz)) prod(sz)
	else              return(NULL)
}

nnz <- function(x) {
    if (length(x) == 0) return(0)
    stopifnot(is.numeric(x) || is.complex(x))
    sum(x != 0)
}

ndims <- function(x) {
	sz <- size(x)
	if (!is.null(sz)) length(sz)
	else              return(NULL)
}

isempty <- function(x) {
	length(x) == 0
}