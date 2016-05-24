def <- function (value, default = "", mode = "character", length.out = NULL)
{
	## Ensure we got a value of a given mode, and if not, use default
	## If length.out is provided, make sure that the returned vector has
	## that length (if needed, cut or recycle 'value')

	## If either NULL or something of length == 0 is in 'value', then,
	## return default
	if (!length(value)) value <- default

	## Coerce to mode...
	res <- switch(as.character(mode[1]),
		logical = as.logical(value),
		character = as.character(value),
		numeric = as.numeric(value),
		double = as.double(value),
		integer = as.integer(value),
		single = as.single(value),
		factor = as.factor(value),
		complex = as.complex(value),
		value)	# This is for unrecognized modes!

	## If length.out is provided, make sure the vector has this length
	if (!is.null(length.out)) {
		if (length(length.out) == 0) length.out <- 1 else
			length.out <- round(as.numeric(length.out[1]))
		res <- rep(res, length.out = length.out)
	}
	res
}
