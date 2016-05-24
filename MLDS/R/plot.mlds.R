`plot.mlds` <-
function(x, standard.scale = FALSE, SD.scale = FALSE, ... ) {
#x, object of class mlds
	par(ask = FALSE)
	f <- if (SD.scale) 2 else 1
	if (standard.scale) {
		ll <- length(x$pscale)
		plot(x$stimulus, x$pscale/x$pscale[ll], ...)
		} else
	plot(x$stimulus, f * x$pscale, ... )
	}

`lines.mlds`<-
function(x, standard.scale = FALSE, SD.scale = FALSE, ... ) {
#x, object of class mlds
	f <- if (SD.scale) 2 else 1
	if (standard.scale) {
		ll <- length(x$pscale)
		lines(x$stimulus, x$pscale/x$pscale[ll], ...)
	} else
	lines(x$stimulus, f * x$pscale, ...)
	}

`points.mlds`<-
function(x, standard.scale = FALSE, SD.scale = FALSE, ... ) {
#x, object of class mlds
	f <- if (SD.scale) 2 else 1
	if (standard.scale) {
		ll <- length(x$pscale)
		points(x$stimulus, x$pscale/x$pscale[ll], ...)
	} else
	points(x$stimulus, f * x$pscale, ...)
	}

`plot.mlbs` <- function(x, standard.scale = FALSE, SD.scale = FALSE, ...){
	 par(ask = FALSE)
	f <- if (SD.scale) 2 else 1
    if (standard.scale) {
        ll <- length(x$pscale)
        plot(x$stimulus, x$pscale/x$pscale[ll], ...)
    }
    else plot(x$stimulus, f * x$pscale, ...)
}

`lines.mlbs` <- function (x, standard.scale = FALSE, SD.scale = FALSE, ...) 
{
	f <- if (SD.scale) 2 else 1
    if (standard.scale) {
        ll <- length(x$pscale)
        lines(x$stimulus, x$pscale/x$pscale[ll], ...)
    }
    else lines(x$stimulus, f * x$pscale, ...)
}

`points.mlbs` <- function (x, standard.scale = FALSE, SD.scale = FALSE, ...) 
{
	f <- if (SD.scale) 2 else 1
    if (standard.scale) {
        ll <- length(x$pscale)
        points(x$stimulus, x$pscale/x$pscale[ll], ...)
    }
    else points(x$stimulus, f * x$pscale, ...)
}
