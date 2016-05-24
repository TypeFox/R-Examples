parse.date <- function(obj){
	varname <- names(obj);
	values <- sapply(obj[[1]]$values, parsevector);
	as.POSIXct(strptime(values, format="%Y-%m-%d %H:%M:%S"));
}
