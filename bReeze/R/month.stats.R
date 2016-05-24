month.stats <-
function(mast, set, signal="v.avg", fun=c("mean", "median", "min", "max", "sd"), subset, digits=3, print=TRUE) {
### calculating monthly statistics

	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	if(missing(set)) set <- "all"
	if(is.null(signal)) stop("Please choose signal")
	if(length(signal)>1) stop("Please choose only one signal")
	if(missing(fun) || length(fun)!=1) fun <- "mean"
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	m.stats.l <- NULL
	years <- unique(mast$timestamp[start:end]$year+1900)
	num.sets <- length(mast$sets)
	unit <- NULL
	
	if(set!="all") { # one set
		if(!is.numeric(set)) set <- match(set, names(mast$sets))
		if(is.na(set)) stop("'set' not found") 
		if(set<0 | set>num.sets) stop("'set' not found")
		if(!any(names(mast$sets[[set]]$data)==signal)) stop("'set' does not contain the choosen signal")
		dat <- mast$sets[[set]]$data[,which(names(mast$sets[[set]]$data)==signal)][start:end]
		m.stats.l <- list(month.stats.int(dat, fun, mast$timestamp[start:end], years, digits))
		names(m.stats.l) <- names(mast$sets)[set]
		unit <- attr(mast$sets[[set]]$data[,which(names(mast$sets[[set]]$data)==signal)], "unit")
	} else { # all sets
		set.index <- NULL
		for(s in 1:num.sets) if(any(names(mast$sets[[s]]$data)==signal)) set.index <- append(set.index, s)
		if(is.null(set.index)) stop("Signal not found in any set")
		
		m.stats.l <- list(month.stats.int(mast$sets[[set.index[1]]]$data[,which(names(mast$sets[[set.index[1]]]$data)==signal)][start:end], fun, mast$timestamp[start:end], years, digits))
		unit <- attr(mast$sets[[set.index[1]]]$data[,which(names(mast$sets[[set.index[1]]]$data)==signal)], "unit")
		
		if(length(set.index) > 1) {
			for(s in 2:length(set.index)) {
				m.stats.df <- month.stats.int(mast$sets[[set.index[s]]]$data[,which(names(mast$sets[[set.index[s]]]$data)==signal)][start:end], fun, mast$timestamp[start:end], years, digits)
				m.stats.l[[length(m.stats.l)+1]] <- m.stats.df
			}
		}
		names(m.stats.l) <- names(mast$sets)[set.index]
	}

	attr(m.stats.l, "unit") <- unit
	attr(m.stats.l, "call") <- list(func="month.stats", mast=deparse(substitute(mast)), set=set, signal=signal, fun=fun, subset=subset, digits=digits, print=print)
	
	class(m.stats.l) <- "month.stats"
	if(print) print(m.stats.l)
	invisible(m.stats.l)
}
