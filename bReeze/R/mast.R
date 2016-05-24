mast <-
function(timestamp, ..., loc=NULL, desc=NULL) {
### creating met mast from several datasets

	if(missing(timestamp)) stop("'timestamp' is mandatory")
	if(!is.null(loc)) {
		if(!is.vector(loc)) stop("'location' must be a vector of latitude and longitude")
		if(length(loc)!=2) stop("'location' must be a numeric vector of latitude and longitude")
		if(!is.numeric(loc)) stop("'location' must be a numeric vector of decimal degrees")
		if(loc[1]>90 || loc[1]<(-90) || loc[2]>180 || loc[2]<(-180)) stop("Coordinates in 'location' out of range - please use decimal degrees")
	}
	
	l <- list(...)
	
	# check sets and timestamp
	num.sets <- length(l)
	if(num.sets<1) stop("No data - please add at least one data set")
	for(i in 1:num.sets) if(class(l[[i]])!="set") stop(names(l)[i], " is no set object")
	if(any(class(timestamp)=="POSIXlt")==FALSE) stop("'timestamp' must be given in POSIXlt format - for reformating use timestamp")
	if(length(timestamp)!=length(l[[1]]$data[,1])) stop("Different length of timestamp and sets")
	
	# check units
	var.names <- NULL
	for(i in 1:num.sets) var.names <- append(var.names, names(l[[i]]$data))
	var.names <- unique(var.names)
	units <- data.frame(matrix(NA, nrow=num.sets, ncol=length(var.names)+1))
	names(units) <- c("height", var.names)
	
	for(i in 1:num.sets) {
		units[i,1] <- attr(l[[i]]$height, "unit")
		for(j in 2:dim(units)[2]) {
			if(!is.null(l[[i]]$data[1,names(units)[j]])) if(!is.null(attr(l[[i]]$data[,names(units)[j]], "unit"))) units[i,j] <- attr(l[[i]]$data[,names(units)[j]], "unit")
		}
	}
	
	for(i in 1:dim(units)[2]) {
		if(length(unique(units[,i])[!is.na(unique(units[,i]))])>1) stop("Different units in ", names(units)[i], " - ", paste(unique(units[,i])[!is.na(unique(units[,i]))], collapse=", "))
	}
	
	# sort by heights
	if(length(l)>1) {
		heights <- c()
		for(i in 1:length(l)) heights <- append(heights, l[[i]]$height)
		height.idx <- cbind(heights, seq(1:length(heights)))
		l <- l[height.idx[order(height.idx[,1], decreasing=TRUE),][,2]]
	}
	
	# name sets
	if(is.null(names(l))) for(i in 1:length(l)) names(l)[i] <- paste0("set", i)
	else for(i in 1:length(l)) if(names(l)[i]=="") names(l)[i] <- paste0("set", i)
	
	if(is.null(loc)) {
		if(is.null(desc)) r <- list(timestamp=timestamp, sets=l)
		else r <- list(timestamp=timestamp, description=desc, sets=l)
	} else {
		if(is.null(desc)) r <- list(timestamp=timestamp, location=loc, sets=l)
		else r <- list(timestamp=timestamp, location=loc, description=desc, sets=l)
	}
	class(r) <- "mast"
	
	return(r)
}
