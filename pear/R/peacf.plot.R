`peacf.plot` <-
function(r)
{
	scale <- 0.5
	type <- attr(r, "type")
	if(is.null(type))
		stop("Error: invalid input. Input must be output from peacf or pepacf"
			)
	if(type == "acf") {
		acf <- scale * r$acf
		sdl <- scale * 1.96 * r$benchmark.sd
		ts.title <- r$title
		ylabel <- "acf"
	}
	else if(type == "pacf") {
		acf <- scale * r$pacf
		sdl <- scale * 1.96 * r$acf.out$benchmark.sd
		ts.title <- r$acf.out$title
		ylabel <- "pacf"
	}
	if(is.null(acf)) {
		stop("error: input argument must be generated from peacf() or pepacf()"
			)
	}
	p <- dim(acf)[[1]]
	lag.max <- dim(acf)[[2]]
	plot(c(1:p, rep(0, lag.max)), c(rep(0, p), 1:lag.max), type = "n", xlab
		 = "period", ylab = "")
	for(imonth in 1:p) {
		segments((imonth - sdl), 0, (imonth - sdl), lag.max)
		segments((imonth + sdl), 0, (imonth + sdl), lag.max)
		for(ilag in 1:lag.max) {
			segments(imonth, ilag, (imonth + acf[imonth, ilag]), 
				ilag)
		}
	}
	title(main = ts.title, ylab = ylabel)
}

