`is.gp` <-
function(x) {
	if (inherits(x, "gp"))
		return (TRUE)
	return (FALSE)
}

