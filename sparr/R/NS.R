NS <- function(data, nstar = NULL, scaler = NA){
    if(class(data)=="data.frame"){
		if(ncol(data)!=2) stop("data.frame must have exactly two columns")
		if(nrow(data)<10) warning("less than 10 observations!")
	} else if(class(data)=="list"){
		if(is.null(data$x)||is.null(data$y)) stop("data list must have two components named 'x' and 'y'")
		if(length(data$x)!=length(data$y)) stop("data components 'x' and 'y' are of unequal lengths")
		if(length(data$x)<10) warning("less than 10 observations!")
		data <- as.data.frame(data)
	} else if(class(data)=="matrix"){
		if(ncol(data)!=2) stop("data matrix must have exactly two columns")
		if(nrow(data)<10) warning("less than 10 observations!")
		data <- as.data.frame(data)
	} else if(class(data)=="ppp"){
		if(is.null(data$x)||is.null(data$y)) stop("data ppp.object must have two non-empty components named 'x' and 'y'")
		if(length(data$x)<10) warning("less than 10 observations!")
		data <- data.frame(cbind(data$x,data$y))
	} else {
		stop("data must be an object of type 'data.frame', 'list', 'matrix', or 'ppp'")
	}
		
	if(is.null(nstar)) n <- nrow(data)
	else n <- nstar
	
	if(n<=0) stop("'nstar' must be positive")
	if(!is.na(scaler)){
		if(scaler<=0) stop("'scaler' must be positive")
	}
	
	if(is.na(scaler)) U <- mean(c(IQR(data[,1])/1.34,IQR(data[,2])/1.34))
	else U <- scaler
    	
    hNS <- U*((1/(2*n))^(1/6))
    return(hNS)
}
