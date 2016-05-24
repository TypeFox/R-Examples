"as.data.frame.lo.smooth" <-
function(x, row.names = NULL, optional = FALSE,...)
{
	d <- dim(x)
	nrows <- d[[1.]]
	dn <- dimnames(x)
	row.names <- dn[[1.]]
	value <- list(x)
	if(length(row.names)) {
		row.names <- as.character(row.names)
		if(length(row.names) != nrows)
			stop(paste("supplied", length(row.names), 
				"names for a data frame with", nrows, "rows"))
	}
	else if(optional)
		row.names <- character(nrows)
	else row.names <- as.character(seq(length = nrows))
	if(!optional)
		names(value) <- deparse(substitute(x))[[1.]]
	attr(value, "row.names") <- row.names
	oldClass(value) <- "data.frame"
	value
}
