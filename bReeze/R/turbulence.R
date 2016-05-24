turbulence <-
function(mast, turb.set, dir.set, num.sectors=12, bins=c(5,10,15,20), subset, digits=3, print=TRUE) {
### calculating mean wind speed and turbulence intensity of sectors

	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	num.sets <- length(mast$sets)
	if(!missing(turb.set) && missing(dir.set)) dir.set <- turb.set
	if(missing(turb.set) && !missing(dir.set)) turb.set <- dir.set
	
	if(!is.numeric(num.sectors)) stop("'num.sectors' must be numeric")
	if(num.sectors<=1) stop("There must be at least 2 sectors")
	if(!is.numeric(turb.set)) turb.set <- match(turb.set, names(mast$sets))
	if(is.na(turb.set)) stop("'turb.set' not found")
	if(turb.set<=0 || turb.set>num.sets) stop("'turb.set' not found")
	if(!is.numeric(dir.set)) dir.set <- match(dir.set, names(mast$sets))
	if(is.na(dir.set)) stop("'dir.set' not found")
	if(dir.set<=0 || dir.set>num.sets) stop("'dir.set' not found")
	if(is.null(mast$sets[[turb.set]]$data$turb.int)) stop("'set' does not contain turbulence intensity data")
	if(is.null(mast$sets[[dir.set]]$data$dir.avg)) stop("'set' does not contain wind direction data")
	if(any(bins<0)) stop("'bins' must be NULL or a vector of positives")
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	v <- mast$sets[[turb.set]]$data$v.avg[start:end]
	tu <- mast$sets[[turb.set]]$data$turb.int[start:end]
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
	
	turb.tbl <- matrix(NA, nrow=num.sectors+1, ncol=num.classes+1)
	# indices for valid data
	idx.val <- !is.na(tu) & !is.na(d)
	idx.v <- !is.na(v)
	
	for(s in 1:num.sectors) {
		# index for direction
		low <- sector.edges[s]
		high <- sector.edges[s+1]
		if(low<high) idx.dir <- d>=low & d<high
		else idx.dir <- d>=low | d<high
		
		if(length(tu[idx.val & idx.dir])<3) turb.tbl[s,1] <- NA
		else turb.tbl[s,1] <- mean(tu[idx.val & idx.dir])
		if(!is.null(bins)) {
			for(c in 1:(num.classes-1)) {
				# index for wind class
				idx.class <- v>=bins[c] & v[start:end]<bins[c+1]
				if(length(tu[idx.val & idx.v & idx.dir & idx.class])<3) turb.tbl[s,c+1] <- NA
				else turb.tbl[s,c+1] <- mean(tu[idx.val & idx.v & idx.dir & idx.class])
			}
			if(length(tu[idx.val & idx.v & idx.dir & v>=bins[num.classes]])<3) turb.tbl[s,num.classes+1] <- NA
			else turb.tbl[s,num.classes+1] <- mean(tu[idx.val & idx.v & idx.dir & v>=bins[num.classes]])
		}
	}
	if(length(tu[start:end])<3) turb.tbl[num.sectors+1,1] <- NA
	else turb.tbl[num.sectors+1,1] <- mean(tu, na.rm=TRUE)
	
	if(!is.null(bins)) {
		for(i in 1:(num.classes-1)) {
			# index for wind class
			idx.class <- v>=bins[i] & v<bins[i+1]
			if(length(tu[idx.val & idx.v & idx.class])<3) turb.tbl[num.sectors+1,i+1] <- NA
			else turb.tbl[num.sectors+1,i+1] <- mean(tu[idx.val & idx.v & idx.class], na.rm=TRUE)
		}
		if(length(tu[idx.val & idx.v & v>=bins[num.classes]])<3) turb.tbl[num.sectors+1,num.classes+1] <- NA
		else turb.tbl[num.sectors+1,num.classes+1] <- mean(tu[idx.val & idx.v & v>=bins[num.classes]])
	}
	
	r.names <- c(paste0("s", 1:num.sectors), "all")
	if(num.sectors==4) r.names <- c("n","e","s","w","all")
	if(num.sectors==8) r.names <- c("n","ne","e","se","s","sw","w","nw","all")
	if(num.sectors==12) r.names <- c("n","nne","ene","e","ese","sse","s","ssw","wsw","w","wnw","nnw","all")
	if(num.sectors==16) r.names <- c("n","nne","ne","ene","e","ese","se","sse","s","ssw","sw","wsw","w","wnw","nw","nnw","all")
	turb.tbl <- data.frame(turb.tbl, row.names=r.names)
	c.names <- c("total")
	if(!is.null(bins)) {
		for(i in 1:(num.classes-1)) c.names <- append(c.names, paste(bins[i], bins[i+1], sep="-"))
		c.names <- append(c.names, paste0(">", bins[num.classes]))
	}
	names(turb.tbl) <- c.names
	
	for(i in 1:length(turb.tbl)) turb.tbl[,i][is.nan(turb.tbl[,i]) | is.na(turb.tbl[,i])] <- 0
	if(sum(turb.tbl[,length(turb.tbl)], na.rm=TRUE)==0) turb.tbl[,length(turb.tbl)] <- NULL
	
	attr(turb.tbl, "call") <- list(func="turbulence", mast=deparse(substitute(mast)), turb.set=turb.set, dir.set=dir.set, num.sectors=num.sectors, bins=bins, subset=subset, digits=digits, print=print)
	turb.tbl <- round(turb.tbl, digits)
	
	class(turb.tbl) <- "turbulence"
	if(print) print(turb.tbl)
	invisible(turb.tbl)
}
