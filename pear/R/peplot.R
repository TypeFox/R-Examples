`peplot` <-
function(z, lag = 1, label = FALSE, mfrow = c(2, 2))
{
	if(2 != length(mfrow)) {
		stop("\nError: mfrow must be vector of length 2. \nExamples:\n     mfrow=c(1,1) produces 1 plot at a time\n     mfrow=c(2,2) produces 4 plots at a time."
			)
	}
	ztitle <- attr(z, "title")
	par(mfrow = mfrow)
	plot.count <- apply(matrix(mfrow, ncol = 2, nrow = 1), MARGIN = 1, FUN
		 = "prod")
	z.title <- attr(z, "title")
	if(is.null(z.title))
		z.title <- " "
	z.names <- attr(z, "period.abb")
	if(is.null(z.names)) {
		z.names <- paste("period", unique(cycle(z)))
	}
	p <- attr(z, "tsp")[3]
	start.month <- cycle(z)[1]
	icount <- 0
	for(imonth in 1:p) {
		icount <- icount + 1
		jmonth <- (imonth - lag - 1) %% p + 1
		y <- z[cycle(z) == imonth]
		x <- z[cycle(z) == jmonth]
		u <- 1:length(y)
		if(start.month != 1) {
			if((imonth >= start.month) && (jmonth < start.month)) 
				{
				y <- y[-1]
				u <- u[-1]
			}
			if((imonth < start.month) && (imonth - lag < 1)) {
				x <- x[ - length(x)]
			}
		}
		if((imonth - lag) < 1) {
			ndrop <- 1 + trunc((lag - imonth)/p)
			n <- length(y)
			y <- y[ - (1:ndrop)]
			u <- u[ - (1:ndrop)]
			x <- x[ - ((n - ndrop + 1):n)]
		}
		plot(x, y, xlab = z.names[jmonth], ylab = z.names[imonth])
		if(label) {
			identify(x, y, labels = u)
		}
		if((icount != p) && (icount %% plot.count) == 0) {
			if(!is.null(ztitle)) {
			        mtext(ztitle, side = 3, outer = TRUE, line = -2, cex = 1.2)
				#par(mfrow = c(1, 1))
				#title(main = ztitle)
				#par(mfrow = mfrow)
			}
			cat("\n Press Enter key for next plot")
			junk <- (as.character(readline()))[1]
		}
		else if((icount %% plot.count) == 0) {
			if(!is.null(ztitle)) {
			        mtext(ztitle, side = 3, outer = TRUE, line = -2, cex = 1.2)
				#par(mfrow = c(1, 1))
				#title(main = ztitle)
				#par(mfrow = mfrow)
			}
		}
	}
	par(mfrow = c(1, 1))
}

