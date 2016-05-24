#.addAxis2------------------------------2013-03-13
# Modified from the function `.addAxis` to 
# add axes to any side of the plot. 
# Note: temporary until incorporation and testing.
#--------------------------------------------NB/RH
.addAxis2 <-
function (side=1:2, xlim, ylim, tckLab, tck, tckMinor, ...) 
{
	tckLab <- rep(tckLab, length.out = 4)
	tck <- rep(tck, length.out = 4)
	tckMinor <- rep(tckMinor, length.out = 4)
	mai <- par()$mai
	# Nick: The following calls to par() are very strange to me as 
	# repeated calls to `.addAxis` sequently reduced the size of cex
	#par(cex = par()$cex * 0.8)
	#par(mai = mai)
	#par(mex = par()$cex)
	lim <- list(xlim,ylim,xlim,ylim)
	rotate <- list(0,90,0,90)
	for (i in side) {
		if (((i %in% c(1,3)) && (par()$xaxt != "n")) || ((i %in% c(2,4)) && (par()$yaxt != "n"))) {
			ticks <- pretty(lim[[i]])
			ticksMinor <- pretty(c(0, diff(ticks)[1]))
			ticksMinor <- (sort(rep(ticks, length(ticksMinor))) + 
				rep(ticksMinor, length(ticks)))
			ticks <- ticks[ticks > lim[[i]][1] & ticks < lim[[i]][2]]
			ticksMinor <- ticksMinor[ticksMinor > lim[[i]][1] & 
				ticksMinor < lim[[i]][2]]
			if (!tckLab[i]) {
				tickLabels <- FALSE
			}
			else {
				tickLabels <- as.character(ticks)
			}
			axis(side = i, at = ticks, labels = tickLabels, tck = tck[i], srt = rotate[[i]], ...)
			axis(side = i, at = ticksMinor, labels = FALSE, tck = tckMinor[i], ...)
		}
	}
	invisible(NULL)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~.addAxis2