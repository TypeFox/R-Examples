triaxes <-
function(lcol = "darkgrey", lty = 2)
{
	lx = 2 / sqrt(3)
	lines(c(lx / 4, lx), c(1 / 2, 0),
	      lty = lty, col = lcol)
	lines(c(3 * lx / 4, 0), c(1 / 2, 0),
		  lty = lty, col = lcol)
	lines(c(lx / 2, lx / 2), c(0, 1),
		  lty = lty, col = lcol)
}

