timestamp <-
function(timestamp, pattern, tz) {
### formatting time stamp (lookup or with given pattern)
	
	if(anyDuplicated(timestamp)) if(any(duplicated(timestamp)==TRUE)) stop("'timestamp' contains duplicates") # sometimes anyDuplicated() founds duplicates although there are no duplicates
	ts <- nts <- NULL
	
	if(missing(tz)) tz <- ""
	if(tz=="?") {
		tz <- ""
		if(grepl("\ [A-Z]+$", timestamp[1])) {
			tz <- tail(strsplit(as.character(timestamp[1]), " ")[[1]], 1)
			message("Time zone found: ", tz)
		}
		else warning("Time zone not recognized - using 'tz=\"\"'", call.=FALSE)
	}
		
	if(missing(pattern)) { # search for pattern
		pattern.list <- read.table(system.file(package="bReeze", "ts_patterns", "patterns.txt"), sep=",")
		pattern <- as.vector(unlist(pattern.list))
				
		for(i in 1:length(pattern)) {
			nts <- strptime(timestamp[1], pattern[i], tz)
			if(is.na(nts)) {
				if(substr(pattern[i], nchar(pattern[i]), nchar(pattern[i]))=="S") nts <- strptime(paste(timestamp[1], "00:00:00"), pattern[i], tz)
				if(substr(pattern[i], nchar(pattern[i]), nchar(pattern[i]))=="M") nts <- strptime(paste(timestamp[1], "00:00"), pattern[i], tz)
			}
			if(!is.na(nts) && substr(nts,1,2)!="00") {
				nts <- strptime(timestamp, pattern[i], tz)
				if(any(is.na(nts)==TRUE)) {
					if(substr(pattern[i], nchar(pattern[i]), nchar(pattern[i]))=="S") nts[which(is.na(nts)==TRUE)] <- strptime(paste(timestamp[which(is.na(nts)==TRUE)], "00:00:00"), pattern[i], tz)
					if(substr(pattern[i], nchar(pattern[i]), nchar(pattern[i]))=="M") nts[which(is.na(nts)==TRUE)] <- strptime(paste(timestamp[which(is.na(nts)==TRUE)], "00:00"), pattern[i], tz)
				}
				if(!any(is.na(nts)==TRUE)) {
					message("Pattern found: ", pattern[i])
					break
				}
			}
		}
		
		if(length(nts)==1 && length(timestamp)!=1) stop("No pattern found")
	} else { # pattern specified
		nts <- strptime(timestamp, pattern, tz)
		if(any(is.na(nts)==TRUE)) {
			if(substr(pattern, nchar(pattern), nchar(pattern))=="S") nts[which(is.na(nts)==TRUE)] <- strptime(paste(timestamp[which(is.na(nts)==TRUE)], "00:00:00"), pattern, tz)
			if(substr(pattern, nchar(pattern), nchar(pattern))=="M") nts[which(is.na(nts)==TRUE)] <- strptime(paste(timestamp[which(is.na(nts)==TRUE)], "00:00"), pattern, tz)
		}
		if(any(is.na(nts)==TRUE) || substr(nts[1],1,2)=="00") stop("Pattern does not match")
	}
	
	return(nts)
}
