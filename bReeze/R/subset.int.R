subset.int <-
function(timestamp, subset) {
###	internal function for data subsets
	
	num.samples <- length(timestamp)
	tz <- attr(timestamp, "tzone")
	if(is.null(tz)) tz <- ""
	if((!any(is.character(subset)) && !any(is.na(subset))) || length(subset)!=2) stop("Please specify 'subset' as vector of start and end time stamp")
	if(is.na(subset[1])) subset[1] <- as.character(timestamp[1])
	if(is.na(subset[2])) subset[2] <- as.character(timestamp[num.samples])
	if(nchar(subset[1])==10) subset[1] <- paste(subset[1], "00:00:00")
	if(nchar(subset[2])==10) subset[2] <- paste(subset[2], "00:00:00")
	start <- strptime(subset[1], "%Y-%m-%d %H:%M:%S", tz[1])
	end <- strptime(subset[2], "%Y-%m-%d %H:%M:%S", tz[1])
	if(is.na(start)) stop("'start' time stamp in 'subset' not correctly formated")
	if(is.na(end)) stop("'end' time stamp in 'subset' not correctly formated")
	if(start<timestamp[1] || start>timestamp[num.samples]) stop("'start' time stamp in 'subset' not in period")
	if(end<timestamp[1] || end>timestamp[num.samples]) stop("'end' time stamp in 'subset' not in period")
	
	match.date <- difftime(timestamp, ISOdatetime(1,1,1,0,0,0), tz=tz[1], units="days") - difftime(start, ISOdatetime(1,1,1,0,0,0), tz=tz[1], units="days")
	start <- which(abs(as.numeric(match.date)) == min(abs(as.numeric(match.date))))	
	match.date <- difftime(timestamp, ISOdatetime(1,1,1,0,0,0), tz=tz[1], units="days") - difftime(end, ISOdatetime(1,1,1,0,0,0), tz=tz[1], units="days")
	end <- which(abs(as.numeric(match.date)) == min(abs(as.numeric(match.date))))	
		
	return(cbind(start, end))
}
