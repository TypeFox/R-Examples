peakabif <- function(abifdata, 
  chanel, 
  npeak, 
  thres = 400/yscale, 
  fig = TRUE,
  chanel.names = c(1:4,105),
  DATA = paste("DATA", chanel.names[chanel], sep = "."),
  tmin = 1/tscale,
  tmax = abifdata$Data[["SCAN.1"]]/tscale,
  tscale = 1000, 
  yscale = 1000,
  irange = (tmin*tscale):(tmax*tscale),
  y = abifdata$Data[[DATA]][irange]/yscale,
  method = "monoH.FC",
  maxrfu = 1000, ...) {
  	
	y[y < thres] <- 0
	heights <- surfaces <- maxis <- starts <- stops <- numeric(npeak)
	
	innoise <- TRUE
	pkidx <- 1
	for (i in 1:length(y)) {
		if (y[i] > 0) {
			if (innoise) {
				starts[pkidx] <- i
				innoise <- FALSE
			}
		 } else {
			if (!innoise) {
				stops[pkidx] <- i - 1
				innoise <- TRUE
				pkidx <- pkidx + 1
			}
		}
	}

	if (fig) par(mfrow = c(4, 4), mar = c(2, 2, 0, 0) + 0.2, oma = c(0,0,2,0))

	for (i in 1:npeak) {
		x <- starts[i]:stops[i]
		if(length(x) <= 2){
			maxis[i] <- NA
			warning("Not all requested peaks were assigned")
			next
		}
		spfun <- splinefun(x, y[x], method = method)
		maxis[i] <- optimize(spfun, interval = range(x), maximum = TRUE)$maximum
		heights[i] <- spfun(maxis[i])
		surfaces[i] <- integrate(spfun, starts[i], stops[i])$value
		if (fig) {
			xx <- (x-1)/tscale + tmin
			plot(xx, y[x], type = "p", las = 1, ylim = range(y), ...)
			abline(h = thres, col = "red")
			lines(xx, spfun(x), col = "blue")
			abline(v = (maxis[i]-1)/tscale + tmin, col = "grey")
		}
	}
	#
	# Compute baseline:
	#
        baseline <- baselineabif(abifdata$Data[[DATA]][irange], maxrfu = maxrfu)
        baseline <- baseline/yscale
	
	if(fig) mtext(paste(deparse(substitute(abifdata)), ",",
	  DATA, ", tmin =", tmin, ", tmax =", tmax, ", thres =", thres, ", npeak =", npeak, ", yscale = ", yscale), side = 3, outer = TRUE)
	
	invisible(list(maxis = (maxis-1) + tmin*tscale, heights = yscale*heights, surfaces = yscale*surfaces, baseline = baseline))
}
