# Calculate Akaike weights
`Weights` <-
function(x)  UseMethod("Weights")

`Weights.model.selection` <-
function(x) {
	i <- type2col(x, "weight")
	structure(item(x, i) / sum(item(x, i)),	names = row.names(x))
}

`Weights.averaging` <-
function(x) {
	x$msTable[, ncol(x$msTable)]
}

`Weights.data.frame` <-
function(x) {
	if(ncol(x) == 2L && colnames(x)[1L] == "df"	&& is.numeric(x[, 2L]))
		return(Weights(x[, 2L]))
	if(ncol(x) == 1L && is.numeric(x[, 1L]))
		return(Weights(x[, 1L]))
	return(NA)
}

`Weights.numeric` <-
function(x) {
	x <- x - min(x)
	d <- exp(-x / 2)
	d / sum(d)
}

`Weights.default` <-
function(x) {
    cry(, "cannot use \"%s\" as 'x'", class(x)[1L])
}
