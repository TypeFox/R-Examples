#'@method print yplocation
#'@S3method print yplocation
print.yplocation <- function(x,...){

	cat("Yplant - location object (class 'yplocation')\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))

	cat("Latitude =",x$lat,"\n")
	cat("Longitude =",x$long,"\n")
	
	cat("Hemisphere (w/e) =",ifelse(x$long<0,"E","W"),"\n")
	cat("Hemisphere (n/s) =",ifelse(x$lat<0,"S","N"),"\n")
	if(x$tzlongset)cat("Nearest timezone border =",x$tzlong,"\n")

}