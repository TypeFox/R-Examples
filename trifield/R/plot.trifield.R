plot.trifield <-
function(x, contours = TRUE,
                         col = topo.colors(256),
                         lab1 = 'A = 0', lab2 = 'B = 0',
                         lab3 = 'C = 0', tribox = TRUE,
                         axis.lines = TRUE, ...)
{
	lx = 2 / sqrt(3)
	plot(range(x$x), range(x$y), type = 'n', asp = 1,
		 axes = F, xlab = NA, ylab = NA, xpd = NA, ...)
	image(x, col = col, add = TRUE)
	if ( axis.lines ) triaxes()
	if ( contours ) contour(x, add = TRUE)
	if ( !is.na(lab1) )
		text(lx / 6 - 0.06, 1 / 3, lab3, srt = 60)
	if ( !is.na(lab2) )
		text(2 * lx / 3 + 0.06, 2 / 3, lab2, srt = -60)
	if ( !is.na(lab3) )
		text(2 * lx / 3, -0.06, lab1)
	if ( tribox ) tribox(lwd = 3)
}

