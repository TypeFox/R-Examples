shadowtext = function(x, y = NULL, labels, col = "white", bg = "grey25", theta = (1:8/4)*pi, r1 = 0.06, r2 = 0.04, ...){
    # produces text with a shadowed border
	xy = xy.coords(x,y)
	xo = r1*strwidth("N")
	yo = r2*strheight("N")
	for (i in theta){
		text(xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ...)
	}
	text(xy$x, xy$y, labels, col=col, ...)
	text(xy$x, xy$y, labels, col=col, ...)
}

