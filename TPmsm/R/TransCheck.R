dataCheck <- function(time1, event1, Stime, event, names, arg=TRUE) {
	if (arg) {Singular <- "Argument"; Plural <- "Arguments";}
	else {Singular <- "Variable"; Plural <- "Variables";}
	Message <- paste(Singular, " '", names[1], "' is not numeric", sep="")
	if ( !is.numeric(time1) ) return(Message)
	Message <- paste(Singular, " '", names[2], "' must be logical or numeric", sep="")
	if ( !( is.logical(event1) | is.numeric(event1) ) ) return(Message)
	Message <- paste(Singular, " '", names[3], "' is not numeric", sep="")
	if ( !is.numeric(Stime) ) return(Message)
	Message <- paste(Singular, " '", names[4], "' must be logical or numeric", sep="")
	if ( !( is.logical(event) | is.numeric(event) ) ) return(Message)
	Message <- paste(Plural, " '", names[1], "', '", names[2], "', '", names[3], " and '", names[4], "' must have the same length", sep="")
	if ( length(time1) != length(event1) | length(time1) != length(Stime) | length(time1) != length(event) ) return(Message)
	Message <- paste(Singular, " '", names[2], "' must be 0 or 1 if numeric and TRUE or FALSE if logical", sep="")
	if ( any( (event1 != 0 & event1 != 1) | (event1 != FALSE & event1 != TRUE) ) ) return(Message)
	Message <- paste(Singular, " '", names[4], "' must be 0 or 1 if numeric and TRUE or FALSE if logical", sep="")
	if ( any( (event != 0 & event != 1) | (event != FALSE & event != TRUE) ) ) return(Message)
	Message <- paste(Plural, " '", names[1], "' and '", names[3], "' must be greater or equal than 0", sep="")
	if ( any(time1 < 0 | Stime < 0) ) return(Message)
	Message <- paste(Singular, " '", names[3], "' must be greater or equal to ", tolower(Singular), " '", names[1], "'", sep="")
	if ( any(Stime < time1) ) return(Message)
	Message <- paste(Plural, " '", names[3], "' and '", names[1], "' must be equal when ", tolower(Singular), " '", names[2], "' equals 0 or FALSE", sep="")
	if ( any(!event1 & Stime != time1) ) return(Message)
	Message <- paste(Singular, " '", names[4], "' must be equal to 0 or FALSE when ", tolower(Singular), " '", names[2], "' equals 0 or FALSE", sep="")
	if ( any(!event1 & event) ) return(Message)
	Message <- paste("When ", tolower(Plural), " '", names[3], "' and '", names[1], "' are equal and ", tolower(Singular), " '", names[2], "' equals 1 or TRUE, ", tolower(Singular), " '", names[4], "' must equal 1 or TRUE", sep="")
	if ( any(time1 == Stime & event1 & !event) ) return(Message)
	return(NULL)
}

TPCheck <- function(object, s, t) {
	if ( !is.survTP(object) ) return("Argument 'object' must be of class 'survTP'")
	if ( !is.numeric(s) ) return("Argument 's' is not numeric")
	if ( !is.numeric(t) ) return("Argument 't' is not numeric")
	if ( !(0 <= s & s <= t) ) return("'s' and 't' must be positive, and s <= t")
	return(NULL)
}

TPCCheck <- function(object, s, t, x) {
	Message <- TPCheck(object, s, t)
	if ( !is.null(Message) ) return(Message)
	if ( !is.numeric(x) ) return("Argument 'x' is not numeric")
	return(NULL)
}

TPWindowCheck <- function(h, nh, ncv, window) {
	if ( !is.numeric(h) ) return("Argument 'h' must be numeric")
	if (length(h) < 1 | length(h) > 4) return("Argument 'h' length must be between 1 and 4")
	if ( any(h <= 0) ) return("Argument 'h' must be greater than 0")
	if ( !( is.numeric(nh) | is.integer(nh) ) ) return("Argument 'nh' must be numeric or integer")
	if (nh <= 1) return("Argument 'nh' must be greater than 1")
	if ( !( is.numeric(ncv) | is.integer(ncv) ) ) return("Argument 'ncv' must be numeric or integer")
	if (ncv < 10) return("Argument 'ncv' must be greater or equal than 10")
	window0 <- c("normal", "epanech", "biweight", "triweight")
	window1 <- c(window0, "box")
	window2 <- c(window1, "tricube", "triangular", "cosine")
	if ( !( window %in% window2 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight', 'triweight', 'box', 'tricube', 'triangular' or 'cosine'")
	return(NULL)
}

TPCWindowCheck <- function(bw, window, method.weights) {
	if ( !( is.character(bw) | is.numeric(bw) ) ) return("Argument 'bw' must be either a character string or a numeric vector")
	if ( is.character(bw) ) {
		if ( !exists(bw, mode="function") ) return( paste("could not find function '", bw, "'", sep="") )
	}
	window0 <- c("normal", "epanech", "biweight", "triweight")
	window1 <- c(window0, "box")
	window2 <- c(window1, "tricube", "triangular", "cosine")
	if (bw %in% c("ALbw", "CVbw", "PBbw") & !( window %in% window0 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight' or 'triweight'")
	else if (bw == "dpik" & !( window %in% window1 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight', 'triweight' or 'box'")
	if ( !( window %in% window2 ) ) return("Argument 'window' must be one of 'normal', 'epanech', 'biweight', 'triweight', 'box', 'tricube', 'triangular' or 'cosine'")
	if ( !method.weights %in% c("NW", "LL") ) return("Argument 'weights' must be one of 'NW' or 'LL'")
	return(NULL)
}

StateCheck <- function(state.names) {
	if (length(state.names) != 3) return("Argument 'state.names' length must be equal to 3")
	if ( length(state.names) != length( unique(state.names) ) ) return("Argument 'state.names' must be unique")
	return(NULL)
}

BootCheck <- function(conf, n.boot, conf.level, method.boot) {
	if ( !is.logical(conf) ) return("Argument 'conf' must be logical")
	if ( !( is.numeric(n.boot) | is.integer(n.boot) ) ) return("Argument 'n.boot' must be numeric or integer")
	if (n.boot <= 1) return("Argument 'n.boot' must be greater than 1")
	if ( !is.numeric(conf.level) ) return("Argument 'conf.level' is not numeric")
	if (conf.level < 0 | conf.level > 1) return("Argument 'conf.level' must be between 0 and 1")
	if ( !( method.boot %in% c("percentile", "basic") ) ) return("Argument 'method.boot' must be one of 'percentile' or 'basic'")
	return(NULL)
}

CvalCheck <- function(boot.cv, cv.full) {
	if ( !is.logical(boot.cv) ) return("Argument 'boot.cv' must be logical")
	if ( !is.logical(cv.full) ) return("Argument 'cv.full' must be logical")
	return(NULL)
}

TPStateBootCheck <- function(object, s, t, state.names, conf, n.boot, conf.level, method.boot) {
	Message <- TPCheck(object, s, t)
	if ( !is.null(Message) ) return(Message)
	Message <- StateCheck(state.names)
	if ( !is.null(Message) ) return(Message)
	Message <- BootCheck(conf, n.boot, conf.level, method.boot)
	return(Message)
}

TPWindowStateBootCvalCheck <- function(object, s, t, h, nh, ncv, window, state.names, conf, n.boot, conf.level, method.boot, boot.cv, cv.full) {
	Message <- TPCheck(object, s, t)
	if ( !is.null(Message) ) return(Message)
	Message <- TPWindowCheck(h, nh, ncv, window)
	if ( !is.null(Message) ) return(Message)
	Message <- StateCheck(state.names)
	if ( !is.null(Message) ) return(Message)
	Message <- BootCheck(conf, n.boot, conf.level, method.boot)
	if ( !is.null(Message) ) return(Message)
	Message <- CvalCheck(boot.cv, cv.full)
	return(Message)
}

TPCWindowStateBootCheck <- function(object, s, t, x, bw, window, method.weights, state.names, conf, n.boot, conf.level, method.boot) {
	Message <- TPCCheck(object, s, t, x)
	if ( !is.null(Message) ) return(Message)
	Message <- TPCWindowCheck(bw, window, method.weights)
	if ( !is.null(Message) ) return(Message)
	Message <- StateCheck(state.names)
	if ( !is.null(Message) ) return(Message)
	Message <- BootCheck(conf, n.boot, conf.level, method.boot)
	return(Message)
}
