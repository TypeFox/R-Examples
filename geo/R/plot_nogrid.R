#' Grid control ?
#' 
#' Grid or border control of some sort?
#' 
#' 
#' @param o o, some outer boundary?
#' @param xgr xgr ?
#' @param ygr ygr ?
#' @param col Color
#' @return No value, some lines added to current geoplot.
#' @note Needs elaboration.
#' @seealso Called by \code{\link{gridaxes}}.
#' @keywords aplot
#' @export plot_nogrid
plot_nogrid <-
function(o, xgr, ygr, col)
{
	frame <- list(x = c(o$x[1], o$x[2], o$x[2], o$x[1], o$x[1]), y = c(
		o$y[1], o$y[1], o$y[2], o$y[2], o$y[1]))
	dx <- (o$x[2] - o$x[1])/100
	ly <- length(ygr)
	lx <- length(xgr)
	lengd <- ly * 2 + lx * 2
	o1 <- o$x[1]
	ind <- c(1:ly)
	my <- mx <- matrix(NA, lengd, 3)
	mx[ind, 1] <- o1
	mx[ind, 2] <- o1 + dx
	my[ind, 1] <- my[ind, 2] <- ygr
	o1 <- o$x[2]
	ind <- c((ly + 1):(ly * 2))
	mx[ind, 1] <- o1 - dx
	mx[ind, 2] <- o1
	my[ind, 1] <- my[ind, 2] <- ygr
	o1 <- o$y[1]
	ind <- c((ly * 2 + 1):(ly * 2 + lx))
	my[ind, 1] <- o1
	my[ind, 2] <- o1 + dx
	mx[ind, 1] <- mx[ind, 2] <- xgr
	o1 <- o$y[2]
	ind <- c((ly * 2 + lx + 1):(ly * 2 + lx * 2))
	my[ind, 1] <- o1 - dx
	my[ind, 2] <- o1
	mx[ind, 1] <- mx[ind, 2] <- xgr
	lines(t(mx), t(my), col = col)
	lines(frame, col = col)
	return(invisible())
}

