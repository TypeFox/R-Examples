error.bars <-
function(x, up, lo, width = 0.02, ...)
{
	xlim <- range(x)
	ci <- diff(xlim) * width
	segments(x, up, x, lo, ...)
	segments(x - ci, up, x + ci, up, ...)
	segments(x - ci, lo, x + ci, lo, ...)
	range(up, lo)
}

