"[.smooth" <-
function(x, ..., drop = FALSE)
{
	cl <- oldClass(x)
	oldClass(x) <- NULL
	ats <- attributes(x)
	ats$dimnames <- NULL
	ats$dim <- NULL
	ats$names <- NULL
	y <- x[..., drop = drop]
	if(!is.null(nas <- ats$NAs)) {
		if(is.null(d <- dim(x)))
			d <- c(length(x), 1.)
		navec <- array(logical(d[1.]), d)
		navec[nas,  ] <- TRUE
		navec <- navec[...]
		nas <- if(is.null(dim(navec))) navec else navec[, 1.]
		nas <- seq(nas)[nas]
		ats$NAs <- nas
	}
	attributes(y) <- c(attributes(y), ats)
	oldClass(y) <- cl
	y
}
