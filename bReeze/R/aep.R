aep <- 
function(profile, pc, hub.h, rho=1.225, avail=1, bins=c(5,10,15,20), sectoral=FALSE, digits=c(3,0,0,3), print=TRUE) {
###	calculating annual energy production
	
	if(missing(profile)) stop("Profile 'profile' is mandatory")
	if(missing(pc)) stop("Power curve 'pc' is mandatory")
	if(missing(hub.h)) stop("Hub heigth 'hub.h' is mandatory")
	if(missing(rho)) rho <- 1.225
	if(missing(avail)) avail <- 1
	if(missing(bins)) bins <- c(5,10,15,20)
	if(missing(sectoral)) sectoral <- FALSE
	if(missing(digits)) digits <- c(3,0,0,3)
	if(missing(print)) print <- TRUE
	
	if(class(profile)!="profile") stop(substitute(profile), " is no profile object")
	if(class(pc)!="pc") stop(substitute(pc), " is no power curve object")
	if(!is.numeric(hub.h)) stop("'hub.h' must be numeric")
	if(!is.numeric(rho)) stop("'rho' must be numeric")
	if(!is.numeric(avail)) stop("'avail' must be numeric")
	if(avail<0 || avail>1) stop("'avail' must be a numeric value between 0 and 1")
	if(any(bins<0)) stop("'bins' must be NULL or a vector of positives")
	if(length(digits)!=4) { 
		digits <- rep(digits, 4)
		warning("'digits' shall be a vector of four values (for wind speed, operation, aep and capacity)", call.=FALSE)
	}
	
	if(is.null(attr(profile, "call")$mast)) stop("Source mast object of ", substitute(profile), " could not be found")
	mast <- get(attr(profile, "call")$mast)
	v.set <- attr(profile, "call")$v.set[1]
	dir.set <- attr(profile, "call")$dir.set
	num.sectors <- attr(profile, "call")$num.sectors
	subset <- attr(profile, "call")$subset
	rho.pc <- attr(pc, "rho")
	rated.p <- attr(pc, "rated.power")
	
	# subset
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	sector.width <- 360/num.sectors
	sectors <- seq(0, 360-sector.width, by=sector.width)
	sector.edges <- c(sectors-sector.width/2, tail(sectors, n=1)+sector.width/2)%%360
	
	v.ref <- mast$sets[[v.set]]$data$v.avg[start:end]
	h.ref <- profile$h.ref 
	dir <- mast$sets[[dir.set]]$data$dir.avg[start:end]
	idx <- !is.na(v.ref) & !is.na(dir)
	v.hh <- v.ref[idx]
	dir <- dir[idx]
	
	if(sectoral) {
		for(i in 1:num.sectors) {
			low <- sector.edges[i]
			high <- sector.edges[i+1]
			if(low<high) sector.idx <- dir>=low & dir<high
			else sector.idx <- dir>=low | dir<high
			v.hh[sector.idx] <- v.hh[sector.idx] * exp(profile$profile[i,"alpha"] * log(hub.h/h.ref))
		}
	} else v.hh <- v.hh * exp(profile$profile["all","alpha"] * log(hub.h/h.ref))
	
	v.max <- max(v.hh, na.rm=TRUE)
	
	if(!is.null(bins)) if(head(bins, 1)!=0) bins <- c(0, bins)
	num.classes <- length(bins)
	if(num.classes>2) {
		for(i in (num.classes-1):2) {
			if(bins[i+1]>=v.max & bins[i]>=v.max) {
				bins <- head(bins, -1)
				num.classes <- length(bins)
			}
		}
	}
	if(!is.null(bins)) if(num.classes==2 && bins[num.classes]>=v.max) stop("Only one wind class found")
	
	aep.tbl <- data.frame(matrix(NA, nrow=num.sectors+1, ncol=num.classes+3))
	r.names <- c(paste0("s", 1:num.sectors),"all")
	if(num.sectors==4) r.names <- c("n","e","s","w","total")
	if(num.sectors==8) r.names <- c("n","ne","e","se","s","sw","w","nw","total")
	if(num.sectors==12) r.names <- c("n","nne","ene","e","ese","sse","s","ssw","wsw","w","wnw","nnw","total")
	if(num.sectors==16) r.names <- c("n","nne","ne","ene","e","ese","se","sse","s","ssw","sw","wsw","w","wnw","nw","nnw","total")
	row.names(aep.tbl) <- r.names
	c.names <- c("wind.speed", "operation", "total")
	if(!is.null(bins)) {
		for(i in 1:(num.classes-1)) c.names <- append(c.names, paste(bins[i], bins[i+1], sep="-"))
		c.names <- append(c.names, paste0(">", bins[num.classes]))
	}
	names(aep.tbl) <- c.names
	
	if(!is.null(bins)) lim.max <- max(5*(trunc(ceiling(max(v.hh, na.rm=TRUE))/5)+1), 5*(trunc(ceiling(max(bins))/5)+1))
	else lim.max <- 5*(trunc(ceiling(max(v.hh, na.rm=TRUE))/5)+1)
		
	for(i in 1:num.sectors) {
		low <- sector.edges[i]
		high <- sector.edges[i+1]
		if(low<high) sector.idx <- dir>=low & dir<high
		else sector.idx <- dir>=low | dir<high
			
		aep.tbl$wind.speed[i] <- round(mean(v.hh[sector.idx], na.rm=TRUE), digits[1])
		aep.tbl$operation[i] <- op <- round(length(v.hh[sector.idx])/length(v.hh)*8760, digits[2])
		wb.par <- weibull.int(v.hh[sector.idx], FALSE)
		if(!is.null(bins)) {
			for(j in 2:num.classes) aep.tbl[i,j+2] <- round(aep.int(wb.par, c(bins[j-1],bins[j]), pc, rho.pc, op, rho, avail), digits[3])
			aep.tbl[i,num.classes+3] <- round(aep.int(wb.par, c(bins[num.classes],lim.max), pc, rho.pc, op, rho, avail), digits[3])
		}
		aep.tbl$total[i] <- round(aep.int(wb.par, c(0,lim.max), pc, rho.pc, op, rho, avail), digits[3])
	}
	aep.tbl$wind.speed[num.sectors+1] <- round(mean(v.hh, na.rm=TRUE), digits=digits[1])
	aep.tbl$operation[num.sectors+1] <- 8760
	
	for(i in 3:(num.classes+3)) aep.tbl[num.sectors+1,i] <- sum(aep.tbl[1:num.sectors,i], na.rm=TRUE)
	for(i in 1:length(aep.tbl)) aep.tbl[,i][is.nan(aep.tbl[,i]) | is.na(aep.tbl[,i])] <- 0
	
	if(!is.null(bins)) if(tail(bins,1)>=v.max) aep.tbl[,length(aep.tbl)] <- NULL
	if(sum(aep.tbl[,length(aep.tbl)], na.rm=TRUE)==0) aep.tbl[,length(aep.tbl)] <- NULL
		
	attr(aep.tbl$wind.speed, "unit") <- "m/s"
	attr(aep.tbl$operation, "unit") <- "h/a"
	attr(aep.tbl$total, "unit") <- "MWh/a"
	
	cap <- round(aep.tbl$total[num.sectors+1] / (rated.p*0.001*8760), digits=digits[4])
	aep <- list(aep=aep.tbl, capacity=cap)
	
	attr(aep, "call") <- list(func="aep", profile=deparse(substitute(profile)), pc=deparse(substitute(pc)), hub.h=hub.h, rho=rho, avail=avail, bins=bins, sectoral=sectoral, digits=digits, print=print)
	class(aep) <- "aep"
	
	if(print) print(aep)
	invisible(aep)
}
