survTP <- function(time1, event1, Stime, event, ...) {
	if ( missing(time1) ) stop("Argument 'time1' is missing, with no default")
	if ( missing(event1) ) stop("Argument 'event1' is missing, with no default")
	if ( missing(Stime) ) stop("Argument 'Stime' is missing, with no default")
	if ( missing(event) ) stop("Argument 'event' is missing, with no default")
	Message <- dataCheck(time1, event1, Stime, event, names=c("time1", "event1", "Stime", "event"), arg=TRUE)
	if ( !is.null(Message) ) stop(Message)
	data <- list("time1"=as.double(time1), "event1"=as.integer(event1), "Stime"=as.double(Stime), "event"=as.integer(event), ...)
	datalen <- length(data)
	if (datalen > 4) {
		datanames <- names(data)
		for (i in 5:datalen) {
			if ( !is.numeric(data[[i]]) ) stop("All additional arguments must be numeric")
			if ( length(data[[i]]) != length(time1) ) stop("All additional arguments must have the same length as arguments 'time1', 'event1', 'Stime' and 'event'")
			if (datanames[i] == "") datanames[i] <- paste("covariate", i-4, sep=".")
			if ( !is.double(data[[i]]) ) data[[i]] <- as.double(data[[i]])
		}
		names(data) <- datanames
	}
	attr(data, "row.names") <- as.integer( 1:length(time1) )
	attr(data, "class") <- "data.frame"
	object <- vector(mode="list", length=1)
	object[[1]] <- data
	attr(object, "class") <- "survTP"
	return(object)
}

is.survTP <- function(x) {
	ret <- inherits(x, "survTP") & is.list(x) & (length(x) >= 1)
	if (!ret) return(ret)
 	ret <- ret & is.data.frame(x[[1]]) & (length(x[[1]]) >= 4)
	if (!ret) return(ret)
	ret <- ret & is.double(x[[1]][[1]]) & is.integer(x[[1]][[2]])
	ret <- ret & is.double(x[[1]][[3]]) & is.integer(x[[1]][[4]])
	for (i in 2:4) ret <- ret & ( length(x[[1]][[1]]) == length(x[[1]][[i]]) )
	if (length(x[[1]]) > 4) {
		for ( i in 5:length(x[[1]]) ) {
			ret <- ret & is.double(x[[1]][[i]])
			ret <- ret & ( length(x[[1]][[1]]) == length(x[[1]][[i]]) )
		}
	}
	ret <- ret & is.integer( attr(x[[1]], "row.names") )
	return(ret)
}
