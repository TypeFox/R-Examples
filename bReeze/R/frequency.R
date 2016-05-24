frequency <-
function(mast, v.set, dir.set, num.sectors=12, bins=c(5,10,15,20), subset, digits=3, print=TRUE) {
### calculating mean wind speed and frequency of sectors

	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	num.sets <- length(mast$sets)
	if(!missing(v.set) && missing(dir.set)) dir.set <- v.set
	if(missing(v.set) && !missing(dir.set)) v.set <- dir.set
	
	if(!is.numeric(v.set)) v.set <- match(v.set, names(mast$sets))
	if(is.na(v.set)) stop("'v.set' not found")
	if(!is.numeric(dir.set)) dir.set <- match(dir.set, names(mast$sets))
	if(is.na(dir.set)) stop("'dir.set' not found")
	
	if(!is.numeric(num.sectors)) stop("'num.sectors' must be numeric")
	if(num.sectors<=1) stop("There must be at least 2 sectors")
	if(v.set<=0 || v.set>num.sets) stop("'v.set' not found")
	if(dir.set<=0 || dir.set>num.sets) stop("'dir.set' not found")
	if(is.null(mast$sets[[v.set]]$data$v.avg)) stop("'set' does not contain average wind speed data")
	if(is.null(mast$sets[[dir.set]]$data$dir.avg)) stop("'set' does not contain wind direction data")
	if(any(bins<0)) stop("'bins' must be NULL or a vector of positives")
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	v <- mast$sets[[v.set]]$data$v.avg[start:end]
	d <- mast$sets[[dir.set]]$data$dir.avg[start:end]
	
	sector.width <- 360/num.sectors
	sectors <- seq(0, 360-sector.width, by=sector.width)
	sector.edges <- c(sectors-sector.width/2, tail(sectors, n=1)+sector.width/2)%%360
	if(!is.null(bins)) if(head(bins, 1)!=0) bins <- c(0, bins)
	num.classes <- length(bins)
	v.max <- max(v, na.rm=TRUE)
	if(num.classes>2) {
		for(i in (num.classes-1):2) {
			if(bins[i+1]>=v.max & bins[i]>=v.max) {
				bins <- head(bins, -1)
				num.classes <- length(bins)
			}
		}
	}
	if(!is.null(bins)) if(num.classes==2 && bins[num.classes]>=v.max) stop("Only one wind class found")
		
	freq.tbl <- matrix(NA, nrow=num.sectors+1, ncol=num.classes+2)
	
	# index for valid data
	idx.val <- !is.na(v) & !is.na(d) & v>=0
	
	for(s in 1:num.sectors) {
		# index for direction
		low <- sector.edges[s]
		high <- sector.edges[s+1]
		if(low<high) idx.dir <- d>=low & d<high
		else idx.dir <- d>=low | d<high
		
		freq.tbl[s,1] <- mean(v[idx.val & idx.dir], na.rm=TRUE)
		freq.tbl[s,2] <- length(v[idx.val & idx.dir]) * 100 / length(d[idx.val])
		if(!is.null(bins)) {
			for(c in 1:(num.classes-1)) {
				# index for wind class
				idx.class <- v>=bins[c] & v<bins[c+1]
				freq.tbl[s,c+2] <- length(v[idx.val & idx.dir & idx.class]) * 100 / length(d[idx.val])
			}
		}
		if(!is.null(bins)) {
			freq.tbl[s,num.classes+2] <- length(v[idx.val & idx.dir & v>=bins[num.classes]]) * 100 / length(d[idx.val])
		}
	}
	freq.tbl[num.sectors+1,1] <- mean(v[idx.val], na.rm=TRUE) # idx.val???
	freq.tbl[num.sectors+1,2] <- sum(freq.tbl[1:num.sectors,2], na.rm=TRUE)
	
	if(!is.null(bins)) for(i in 3:(num.classes+2)) freq.tbl[num.sectors+1,i] <- sum(freq.tbl[1:num.sectors,i], na.rm=TRUE)
	
	r.names <- c(paste0("s", 1:num.sectors),"all")
	if(num.sectors==4) r.names <- c("n","e","s","w","all")
	if(num.sectors==8) r.names <- c("n","ne","e","se","s","sw","w","nw","all")
	if(num.sectors==12) r.names <- c("n","nne","ene","e","ese","sse","s","ssw","wsw","w","wnw","nnw","all")
	if(num.sectors==16) r.names <- c("n","nne","ne","ene","e","ese","se","sse","s","ssw","sw","wsw","w","wnw","nw","nnw","all")
	freq.tbl <- data.frame(freq.tbl, row.names=r.names)
	c.names <- c("wind.speed","total")
	if(!is.null(bins)) {
		for(i in 1:(num.classes-1)) c.names <- append(c.names, paste(bins[i], bins[i+1], sep="-"))
		c.names <- append(c.names, paste0(">", bins[num.classes]))
	}
	names(freq.tbl) <- c.names
	
	for(i in 1:length(freq.tbl)) freq.tbl[,i][is.nan(freq.tbl[,i]) | is.na(freq.tbl[,i])] <- 0
	if(sum(freq.tbl[,length(freq.tbl)], na.rm=TRUE)==0) freq.tbl[,length(freq.tbl)] <- NULL
	
	unit <- attr(freq.tbl, "units") <- c(attr(mast$sets[[v.set]]$data$v.avg, "unit"), "%")
	attr(freq.tbl, "call") <- list(func="frequency", mast=deparse(substitute(mast)), v.set=v.set, dir.set=dir.set, num.sectors=num.sectors, bins=bins, subset=subset, digits=digits, print=print)
	freq.tbl <- round(freq.tbl, digits)
	
	class(freq.tbl) <- "frequency"
	if(print) print(freq.tbl)
	invisible(freq.tbl)
}
