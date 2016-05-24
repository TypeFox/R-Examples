tribox <-
function(...)
{
	lx = 2 / sqrt(3)
	tribox = cbind(c(0, lx / 2, lx), c(0, 1, 0))
	polygon(tribox, ...)
}

