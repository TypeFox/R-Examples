availability <- 
function(mast, v.set, dir.set, subset, digits=1, print=TRUE) {
### check availability for pairs of windspeed and direction - effective data period
	
	if(class(mast)!="mast") stop(substitute(mast), " is no mast object")
	num.sets <- length(mast$sets)
	if(missing(v.set) && missing(dir.set)) v.set <- "all"
	if(!missing(v.set) && missing(dir.set)) dir.set <- v.set
	if(missing(v.set) && !missing(dir.set)) v.set <- dir.set
	
	if(!is.numeric(v.set)) if(!(length(v.set)==1 && any(v.set=="all"))) v.set <- match(v.set, names(mast$sets))
	if(any(is.na(v.set))) stop("'v.set' not found")
	if(!is.numeric(dir.set)) if(!(length(dir.set)==1 && any(dir.set=="all"))) dir.set <- match(dir.set, names(mast$sets))
	if(any(is.na(dir.set))) stop("'dir.set' not found")
	
	# subset
	if(missing(subset)) subset <- c(NA, NA)
	start.end <- subset.int(mast$timestamp, subset)
	start <- start.end[1]
	end <- start.end[2]
	
	num.samples <- length(mast$timestamp[start:end])
	start.year <- mast$timestamp[start:end]$year[1]+1900
	end.year <- mast$timestamp[start:end]$year[num.samples]+1900
	start.month <- mast$timestamp[start:end]$mon[1]+1
	end.month <- mast$timestamp[start:end]$mon[num.samples]+1
	start.day <- mast$timestamp[start:end]$mday[1]
	end.day <- mast$timestamp[start:end]$mday[num.samples]
	if(start.year==end.year) num.months <- end.month-start.month+1
	if(start.year!=end.year) num.months <- 13-start.month+end.month + 12*(end.year-start.year-1)
	period.days <- as.numeric(mast$timestamp[start:end][num.samples]-mast$timestamp[start:end][1])
	
	if(all(v.set!="all")) {
		if(length(v.set)==1 && length(dir.set)==1) { # one set
			if(v.set<0 || v.set>num.sets) stop("'v.set' not found")
			if(dir.set<0 || dir.set>num.sets) stop("'dir.set' not found")
			if(is.null(mast$sets[[v.set]]$data$v.avg)) stop("'v.set' does not contain average wind speed data")
			if(is.null(mast$sets[[dir.set]]$data$dir.avg)) stop("'dir.set' does not contain average wind direction data")
			if(any(attr(mast$sets[[v.set]]$data, "clean")=="v.avg") && any(attr(mast$sets[[dir.set]]$data, "clean")=="dir.avg")) message("Set(s) not cleaned - cleaning of wind speed v.avg and wind direction dir.avg using 'clean' is recommended to avoid overestimated availability")
			avail <- list(availability.int(mast$sets[[v.set]]$data$v.avg[start:end], mast$sets[[dir.set]]$data$dir.avg[start:end], mast$timestamp[start:end], start.year, start.month, num.months, period.days, digits))
			if(v.set==dir.set) names(avail) <- names(mast$sets)[v.set]
			else names(avail) <- paste0(names(mast$sets)[v.set], "_", names(mast$sets)[dir.set])
		} else { # list of sets
			if(length(v.set)==length(dir.set)) {
				avail <- total <- NULL
				uncleaned <- 0
				
				for(s in 1:length(v.set)) { # x/x
					if(v.set[s]<0 || v.set[s]>num.sets) stop("'v.set' not found")
					if(dir.set[s]<0 || dir.set[s]>num.sets) stop("'dir.set' not found")
					if(is.null(mast$sets[[v.set[s]]]$data$v.avg)) stop("'v.set' does not contain average wind speed data")
					if(is.null(mast$sets[[dir.set[s]]]$data$dir.avg)) stop("'dir.set' does not contain average wind direction data")
					if(any(attr(mast$sets[[v.set[s]]]$data, "clean")=="v.avg") && any(attr(mast$sets[[dir.set[s]]]$data, "clean")=="dir.avg")) message("Set(s) not cleaned - cleaning of wind speed v.avg and wind direction dir.avg using 'clean' is recommended to avoid overestimated availability")
					
					avail.s <- availability.int(mast$sets[[v.set[s]]]$data$v.avg[start:end], mast$sets[[dir.set[s]]]$data$dir.avg[start:end], mast$timestamp[start:end], start.year, start.month, num.months, period.days, digits)
					if(!is.null(avail)) avail[[length(avail)+1]] <- avail.s
					if(is.null(avail)) avail <- list(avail.s)
					if(any(attr(mast$sets[[v.set[s]]]$data, "clean")=="v.avg") && any(attr(mast$sets[[v.set[s]]]$data, "clean")=="dir.avg")) uncleaned <- uncleaned+1
					if(v.set[s]==dir.set[s]) names(avail)[s] <- names(mast$sets)[v.set[s]]
					else names(avail)[s] <- paste0(names(mast$sets)[v.set[s]], "_", names(mast$sets)[dir.set[s]])
				}
			} else if(length(v.set)!=length(dir.set) && (length(v.set)==1 || length(dir.set)==1)) { # x/1 or 1/x
				avail <- total <- NULL
				uncleaned <- 0
				
				if(length(v.set)==1) {
					for(s in 1:length(dir.set)) {
						if(v.set<0 || v.set>num.sets) stop("'v.set' not found")
						if(dir.set[s]<0 || dir.set[s]>num.sets) stop("'dir.set' not found")
						if(is.null(mast$sets[[v.set]]$data$v.avg)) stop("'v.set' does not contain average wind speed data")
						if(is.null(mast$sets[[dir.set[s]]]$data$dir.avg)) stop("'dir.set' does not contain average wind direction data")
						if(any(attr(mast$sets[[v.set]]$data, "clean")=="v.avg") && any(attr(mast$sets[[dir.set[s]]]$data, "clean")=="dir.avg")) message("Set(s) not cleaned - cleaning of wind speed v.avg and wind direction dir.avg using 'clean' is recommended to avoid overestimated availability")
						
						avail.s <- availability.int(mast$sets[[v.set]]$data$v.avg[start:end], mast$sets[[dir.set[s]]]$data$dir.avg[start:end], mast$timestamp[start:end], start.year, start.month, num.months, period.days, digits)
						if(!is.null(avail)) avail[[length(avail)+1]] <- avail.s
						if(is.null(avail)) avail <- list(avail.s)
						if(any(attr(mast$sets[[v.set]]$data, "clean")=="v.avg") && any(attr(mast$sets[[v.set]]$data, "clean")=="dir.avg")) uncleaned <- uncleaned+1
						if(v.set==dir.set[s]) names(avail)[s] <- names(mast$sets)[v.set]
						else names(avail)[s] <- paste0(names(mast$sets)[v.set], "_", names(mast$sets)[dir.set[s]])
					}
				} else {
					for(s in 1:length(v.set)) {
						if(v.set[s]<0 || v.set[s]>num.sets) stop("'v.set' not found")
						if(dir.set<0 || dir.set>num.sets) stop("'dir.set' not found")
						if(is.null(mast$sets[[v.set[s]]]$data$v.avg)) stop("'v.set' does not contain average wind speed data")
						if(is.null(mast$sets[[dir.set]]$data$dir.avg)) stop("'dir.set' does not contain average wind direction data")
						if(any(attr(mast$sets[[v.set[s]]]$data, "clean")=="v.avg") && any(attr(mast$sets[[dir.set]]$data, "clean")=="dir.avg")) message("Set(s) not cleaned - cleaning of wind speed v.avg and wind direction dir.avg using 'clean' is recommended to avoid overestimated availability")
						
						avail.s <- availability.int(mast$sets[[v.set[s]]]$data$v.avg[start:end], mast$sets[[dir.set]]$data$dir.avg[start:end], mast$timestamp[start:end], start.year, start.month, num.months, period.days, digits)
						if(!is.null(avail)) avail[[length(avail)+1]] <- avail.s
						if(is.null(avail)) avail <- list(avail.s)
						if(any(attr(mast$sets[[v.set[s]]]$data, "clean")=="v.avg") && any(attr(mast$sets[[v.set[s]]]$data, "clean")=="dir.avg")) uncleaned <- uncleaned+1
						if(v.set[s]==dir.set) names(avail)[s] <- names(mast$sets)[v.set[s]]
						else names(avail)[s] <- paste0(names(mast$sets)[v.set[s]], "_", names(mast$sets)[dir.set])
					}
				}
			}
		}
	} else { # all sets
		set.index <- NULL
		for(s in 1:num.sets) if(!is.null(mast$sets[[s]]$data$v.avg) && !is.null(mast$sets[[s]]$data$dir.avg)) set.index <- append(set.index, s)
		if(is.null(set.index)) stop("No pairs of wind speed and wind direction data found")
		avail <- total <- NULL
		uncleaned <- 0
		
		for(s in 1:length(set.index)) {
			avail.s <- availability.int(mast$sets[[set.index[s]]]$data$v.avg[start:end], mast$sets[[set.index[s]]]$data$dir.avg[start:end], mast$timestamp[start:end], start.year, start.month, num.months, period.days, digits)
			if(!is.null(avail)) avail[[length(avail)+1]] <- avail.s
			if(is.null(avail)) avail <- list(avail.s)
			if(any(attr(mast$sets[[set.index[s]]]$data, "clean")=="v.avg") && any(attr(mast$sets[[set.index[s]]]$data, "clean")=="dir.avg")) uncleaned <- uncleaned+1
		}
		names(avail) <- names(mast$sets)[set.index]
		if(uncleaned>0) message(uncleaned, " of ", length(set.index), " sets were not cleaned - cleaning of wind speed v.avg and wind direction dir.avg using 'clean' is recommended to avoid overestimated availability")
	}
	
	attr(avail, "call") <- list(func="availability", mast=deparse(substitute(mast)), v.set=v.set, dir.set=dir.set, subset=subset, digits=digits, print=print)
	
	class(avail) <- "availability"
	if(print) print(avail)
	invisible(avail)
}
