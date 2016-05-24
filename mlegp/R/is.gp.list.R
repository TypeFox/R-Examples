`is.gp.list` <-
function(x) {
	if (inherits(x, "gp.list"))
		return (TRUE)
	return (FALSE)
}

