"butterfly" <-
function (loc, color, scale=.1, rand=.1) 
{
	x <- runif(2, -scale, scale)
	y <- runif(2, -scale, scale)
	slope <- diff(y)/diff(x)
	orslope <- -1/slope
	len <- sqrt(diff(x)^2 + diff(y)^2 * runif(1, 1-rand, 1))
	p1x <- x[1] - len/sqrt(1 + orslope^2)
	p1y <- orslope * (p1x - x[1]) + y[1]
	p2x <- x[2] - len/sqrt(1 + orslope^2)
	p2y <- orslope * (p2x - x[2]) + y[2]
	polygon(c(x, p1x, p2x) + loc[1], c(y, p1y, p2y) + loc[2], col=color,
		border=NA)

}

