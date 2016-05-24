plotCorCoeff <- function(m, idx, col=c("black"), breaks="Scott", ...) {
	col.default <- "grey75"
	
	if(length(col) != length(names(idx))) {
		stop("The number of colours (", length(col) ,") does not match the number of indices list elements (", length(names(idx)) ,").")
	}
	
	# get data from just one triangle and only that which is not NA
	dat <- m[upper.tri(m) & !is.na(m)]
	all.hist <- hist(dat, plot=FALSE, breaks=breaks)
	
	# get the break information from all.hist
	breaks <- all.hist$breaks
	
	# plot the distribution for all m values
	plot.new()
	plot.window(xlim=range(breaks),
			ylim=range(0, all.hist$counts))
	rect(breaks[-length(breaks)], 0,
			breaks[-1], all.hist$counts, col=col.default, ...)
	
	# for each vector of indices in the idx list, superimpose a distribution with the same breaks as all.hist
	for (i in 1:length(idx)) {
		dat <- m[intersect(idx[[i]], which(upper.tri(m)))]
		
		i.hist <- hist(dat, plot=FALSE, breaks=breaks)
		rect(i.hist$breaks[-length(i.hist$breaks)], 0,
				i.hist$breaks[-1], i.hist$counts, col=col[i], ...)
		
	}
	
	axis(1)
	axis(2)
	title(main="Density Distribution of Correlation Coefficients", xlab="Correlation Coefficient", ylab="Frequency")
	
	# Still need to plot the legend
	if( length(names(idx)) >= 1 ) {
		legend(-1, max(all.hist$counts), fill=c(col.default, col), legend=c("All", names(idx)))
	}
}
