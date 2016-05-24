iqr <- function(x){ return(diff(quantile(x,c(.25,.75),na.rm=T))) }

# dhist <- function(x, a=10*iqr(x),
dhist <- function(x, a = 5*iqr(x), nbins=nclass.Sturges(x), rx = range(x, na.rm=TRUE), eps=.15, xlab = "x", plot = TRUE, lab.spikes=TRUE) {
#x is the data
#a is the scaling factor, default is 5 * IQR
# nbins is the number of bins, default is assigned by the Stuges method
#rx is the range used for the left of the left-most bin to the right of the
#	right-most bin
#xlab is label for the x axis
#plot = TRUE produces the plot, F returns the heights, breaks and counts
#lab.spikes = TRUE labels the % of data in the spikes

if(is.character(nbins))
        nbins <- switch(casefold(nbins),
                sturges = nclass.Sturges(x),
                fd = nclass.FD(x),
                scott = nclass.scott(x),
                stop("Nclass method not recognized"))
        else if(is.function(nbins))
                nbins <- nbins(x)

	x <- sort(x[!is.na(x)])
	if(a == 0)
		a <- diff(range(x))/100000000
	if(a != 0 & a != Inf) {
		n <- length(x)
		h <- (rx[2] + a - rx[1])/nbins
		ybr <- rx[1] + h * (0:nbins)
		yupper <- x + (a * (1:n))/n
		# upper and lower corners in the ecdf
		ylower <- yupper - a/n
		#
		cmtx <- cbind(cut(yupper, breaks = ybr), cut(yupper, breaks = 
			ybr, left.include = T), cut(ylower, breaks = ybr),
			cut(ylower, breaks = ybr, left.include = T))
		cmtx[1, 3] <- cmtx[1, 4] <- 1
		# to replace NAs when default r is used
		cmtx[n, 1] <- cmtx[n, 2] <- nbins
		#
		#checksum <- apply(cmtx, 1, sum) %% 4
		checksum <- (cmtx[, 1] + cmtx[, 2] + cmtx[, 3] + cmtx[, 4]) %%
			4
		# will be 2 for obs. that straddle two bins
		straddlers <- (1:n)[checksum == 2]
		# to allow for zero counts
		if(length(straddlers) > 0) counts <- table(c(1:nbins, cmtx[
				 - straddlers, 1])) else counts <- table(c(
				1:nbins, cmtx[, 1]))
		counts <- counts - 1
		#
		if(length(straddlers) > 0) {
			for(i in straddlers) {
				binno <- cmtx[i, 1]
				theta <- ((yupper[i] - ybr[binno]) * n)/a
				counts[binno - 1] <- counts[binno - 1] + (
					1 - theta)
				counts[binno] <- counts[binno] + theta
			}
		}
		xbr <- ybr
		xbr[-1] <- ybr[-1] - (a * cumsum(counts))/n
		spike<-eps*diff(rx)/nbins
		flag.vec<-c(diff(xbr)<spike,F)
		if ( sum(abs(diff(xbr))<=spike) >1) {
		xbr.new<-xbr
		counts.new<-counts
		diff.xbr<-abs(diff(xbr))
		amt.spike<-diff.xbr[length(diff.xbr)]
		for (i in rev(2:length(diff.xbr))) {
			if (diff.xbr[i-1]<=spike&diff.xbr[i]<=spike&
				!is.na(diff.xbr[i])) {
				amt.spike<-amt.spike+diff.xbr[i-1]
				counts.new[i-1]<-counts.new[i-1]+counts.new[i]
				xbr.new[i]<-NA
				counts.new[i]<-NA
				flag.vec[i-1]<-T
			}
			else amt.spike<-diff.xbr[i-1]
		}
		flag.vec<-flag.vec[!is.na(xbr.new)]
		flag.vec<-flag.vec[-length(flag.vec)]
		counts<-counts.new[!is.na(counts.new)]
		xbr<-xbr.new[!is.na(xbr.new)]

		}
		else flag.vec<-flag.vec[-length(flag.vec)]
		widths <- abs(diff(xbr))
		# N.B. argument "widths" in barplot must be xbr
		heights <- counts/widths
	}
		bin.size <- length(x)/nbins
		cut.pt <- unique(c(min(x) - abs(min(x))/1000, approx(seq(length(
			x)), x, (1:(nbins - 1)) * bin.size, rule = 2)$y, max(
			x)))
		aa <- hist(x, breaks = cut.pt, plot = FALSE)
		# aa <- hist(x, breaks = cut.pt, plot = FALSE, probability = TRUE)
	if(a == Inf) {
		heights <- aa$counts
		xbr <- aa$breaks
	}
	amt.height<-3
	q75<-quantile(heights,.75)
	if (sum(flag.vec)!=0) {
		amt<-max(heights[!flag.vec])
		ylim.height<-amt*amt.height
		ind.h<-flag.vec&heights> ylim.height
		flag.vec[heights<ylim.height*(amt.height-1)/amt.height]<-F
		heights[ind.h] <- ylim.height
	}
	amt.txt<-0
	end.y<-(-10000)
	if(plot) {
		barplot(heights, abs(diff(xbr)), space = 0, density = -1, xlab = 
			xlab, plot = TRUE, xaxt = "n",yaxt='n')
		at <- pretty(xbr)
		axis(1, at = at - xbr[1], labels = as.character(at))
        if (lab.spikes) {
                if (sum(flag.vec)>=1) {
                usr<-par('usr')
		for ( i in seq(length(xbr)-1)) {
		if (!flag.vec[i]) {
			amt.txt<-0
			if (xbr[i]-xbr[1]<end.y) amt.txt<-1
			}
		else {
			amt.txt<-amt.txt+1
			end.y<-xbr[i]-xbr[1]+3*par('cxy')[1]
			}
			if (flag.vec[i]) {
			txt<-paste(' ',format(round(counts[i]/
				sum(counts)*100)),'%',sep='')
			par(xpd=T)
			text(xbr[i+1]-xbr[1],ylim.height-par('cxy')[2]*(amt.txt-1),txt, adj=0)
}}
                }
                else print('no spikes or more than one spike')
        }
		invisible(list(heights = heights, xbr = xbr))
	}
	else {
		return(list(heights = heights, xbr = xbr,counts=counts))
	}
}
