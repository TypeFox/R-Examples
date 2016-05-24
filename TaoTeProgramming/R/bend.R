"bend" <-
function (num=2e4, xdelta=100, ydelta=200, sd=1)
{
	x <- cumsum(rnorm(num, sd=sd) + cos(seq(0, runif(1) * pi, length=num) 
		* xdelta))
	y <- cumsum(rnorm(num, sd=sd) + sin(seq(0, runif(1) * pi, length=num) 
		* ydelta))
	cbind(x, y)
}

