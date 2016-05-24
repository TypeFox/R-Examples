"batons2" <- 
function(..., waist = FALSE){
	bstats <- boxplot(..., plot=FALSE)
	n <- ncol(bstats$stats)
	max.range <- range(unlist(bstats[c(1,3:4)]))
	# start plot
	cl <- match.call()
	if(is.na(match("xlim", names(cl)))){xl <- c(0, n+1)}
	else{xl <- cl$xlim}
	if(is.na(match("ylim", names(cl)))){yl <- max.range + diff(max.range)*c(-0.1, 0.1)}
	else{yl <- cl$ylim}
	plot.default(x = seq(n), y = max.range, type="n", xlim=xl, ylim=yl, ...)
	# draw extremes
	segments(seq(n), bstats$stats[1,], seq(n), bstats$stats[5,])
}