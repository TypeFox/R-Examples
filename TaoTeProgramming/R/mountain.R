"mountain" <-
function (values, color, df=7) 
{
	x <- seq(0, 1, length=length(values))
	top <- smooth.spline(x, values, df=df)
	polygon(c(top$x, 1, 0), c(pmax(top$y, 0), 0, 0), col=color, border=NA)
}

