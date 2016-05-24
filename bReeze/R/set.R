set <- 
function(height, desc, v.avg, v.max, v.min, v.std, dir.avg, dir.std, tmp, ...) {
### creating dataset in particular height

	if(missing(height)) stop("'height' is mandatory")
	if(!is.numeric(height)) stop("'height' must be numeric")
	if(missing(v.avg) && missing(v.max) && missing(v.min) && missing(v.std) && missing(dir.avg) && missing(dir.std) && missing(...)) stop("No data")
	if(missing(desc)) desc <- NULL
	if(missing(v.avg)) v.avg <- NULL
	if(missing(v.max)) v.max <- NULL
	if(missing(v.min)) v.min <- NULL
	if(missing(v.std)) v.std <- NULL
	if(missing(dir.avg)) dir.avg <- NULL
	if(missing(dir.std)) dir.std <- NULL
	if(missing(tmp)) tmp <- NULL
		
	if(!is.null(v.std) && !is.null(v.avg)) turb.int <- v.std / v.avg
	else turb.int <- NULL
	
	data <- data.frame(cbind(v.avg, v.max, v.min, v.std, dir.avg, dir.std, turb.int, ...))
	attr(height, "unit") <- "m"
	if(!is.null(data$v.avg)) attr(data$v.avg, "unit") <- "m/s"
	if(!is.null(data$v.max)) attr(data$v.max, "unit") <- "m/s"
	if(!is.null(data$v.min)) attr(data$v.min, "unit") <- "m/s"
	if(!is.null(data$v.std)) attr(data$v.std, "unit") <- "m/s"
	if(!is.null(data$dir.avg)) attr(data$dir.avg, "unit") <- "deg"
	if(!is.null(data$dir.std)) attr(data$dir.std, "unit") <- "deg"
	if(!is.null(data$tmp)) attr(data$tmp, "unit") <- "deg C"
	if(!is.null(data$turb.int)) attr(data$turb.int, "unit") <- "-"
	
	if(is.null(desc)) l <- list(height=height, data=data)
	else l <- list(height=height, description=desc, data=data)
	class(l) <- "set"
	
	return(l)
}
