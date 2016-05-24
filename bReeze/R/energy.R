energy <-
function(wb, rho=1.225, bins=c(5,10,15,20), digits=0, print=TRUE) {
###	calculating wind energy per sector
	
	if(class(wb)!="weibull") stop(substitute(wb), " is no weibull object")
	if(any(bins<0)) stop("'bins' must be NULL or a vector of positives")

	if(is.null(attr(wb, "call")$mast)) stop("Source mast object of ", substitute(wb), " could not be found")
	mast <- get(attr(wb, "call")$mast)
	v.set <- attr(wb, "call")$v.set
	dir.set <- attr(wb, "call")$dir.set
	num.sectors <- attr(wb, "call")$num.sectors
	subset <- attr(wb, "call")$subset
	
	# subset
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	lim <- c(0, 5*(trunc(ceiling(max(mast$sets[[v.set]]$data$v.avg[start:end], na.rm=TRUE))/5)+1))
	
	if(!is.null(bins)) if(head(bins, 1)!=0) bins <- c(0, bins)
	num.classes <- length(bins)
	v.max <- max(mast$sets[[v.set]]$data$v.avg[start:end], na.rm=TRUE)
	if(num.classes>2) {
		for(i in (num.classes-1):2) {
			if(bins[i+1]>=v.max && bins[i]>=v.max) {
				bins <- head(bins, -1)
				num.classes <- length(bins)
			}
		}
	}
	if(!is.null(bins)) if(num.classes==2 && bins[num.classes]>=v.max) stop("Only one wind class found")

	energy.tbl <- data.frame(matrix(NA, nrow=num.sectors+1, ncol=num.classes+1))
	r.names <- c(paste0("s", 1:num.sectors),"all")
	if(num.sectors==4) r.names <- c("n","e","s","w","all")
	if(num.sectors==8) r.names <- c("n","ne","e","se","s","sw","w","nw","all")
	if(num.sectors==12) r.names <- c("n","nne","ene","e","ese","sse","s","ssw","wsw","w","wnw","nnw","all")
	if(num.sectors==16) r.names <- c("n","nne","ne","ene","e","ese","se","sse","s","ssw","sw","wsw","w","wnw","nw","nnw","all")
	row.names(energy.tbl) <- r.names
	c.names <- c("total")
	if(!is.null(bins)) {
		for(i in 1:(num.classes-1)) c.names <- append(c.names, paste(bins[i], bins[i+1], sep="-"))
		c.names <- append(c.names, paste0(">", bins[num.classes]))
	}
	names(energy.tbl) <- c.names
	
	freq.l <- frequency(mast, v.set, dir.set, num.sectors, bins, subset, print=FALSE)
	freq <- freq.l[[2]]
	if(length(freq.l)>2) for(i in 3:length(freq.l)) freq <- cbind(freq, freq.l[[i]])
	if(!is.null(bins)) freq <- data.frame(freq)
	
	for(i in 1:num.sectors) {
		if(!is.null(bins)) {
			for(j in 2:dim(freq)[2]) energy.tbl[i,j] <- round(energy.int(lim, wb$k[i], wb$A[i], rho)*freq[i,j]/100, digits)
			energy.tbl[i,1] <- round(energy.int(lim, wb$k[i], wb$A[i], rho)*freq[i,1]/100, digits)
		} else {
			energy.tbl[i,1] <- round(energy.int(lim, wb$k[i], wb$A[i], rho)*freq[i]/100, digits)
		}
	}
	for(i in 1:(num.classes+1)) energy.tbl[num.sectors+1,i] <- sum(energy.tbl[1:num.sectors,i], na.rm=TRUE)
	for(i in 1:length(energy.tbl)) energy.tbl[,i][is.nan(energy.tbl[,i]) | is.na(energy.tbl[,i])] <- 0
	
	if(!is.null(bins)) if(tail(bins,1)>=v.max) energy.tbl[,length(energy.tbl)] <- NULL
	if(sum(energy.tbl[,length(energy.tbl)], na.rm=TRUE)==0) energy.tbl[,length(energy.tbl)] <- NULL
	
	attr(energy.tbl, "unit") <- "kWh/m^2/a"	
	attr(energy.tbl, "call") <- list(func="energy", wb=deparse(substitute(wb)), rho=rho, bins=bins, digits=digits, print=print)
	energy.tbl <- round(energy.tbl, digits)
	
	class(energy.tbl) <- "energy"
	if(print) print(energy.tbl)
	invisible(energy.tbl)
}
