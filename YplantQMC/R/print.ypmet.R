#'@method print ypmet
#'@S3method print ypmet
print.ypmet <- function(x,ndigits=2,...){

	cat("Yplant - daily weather object (class 'ypmet').\n\n")
	cat(paste(c(rep("-",30),"\n"),collapse=""))
	
	mayberound <- function(x,n=2)ifelse(is.character(x),x,round(x,n))
	
	if(x$method == "input")
		cat("Weather data (partially) input by user.\n\n")
	if(!(x$daylength == "not calculated")){
	cat("Month                :", x$month,"\n")
	cat("Day of month         :", x$day,"\n")
	if(x$month == 5 & x$day == 22)cat("***Remko's birthday!***\n")
	cat("Day length           :", mayberound(x$daylength,2),"hours.\n")
	cat("Sun rise             :", mayberound(x$sunrise,2),"hours.\n")
	cat("Sun set              :", mayberound(x$sunset,2),"hours.\n")
	}
	cat("Number of time steps :", nrow(x$dat),"\n\n")
	cat("Data :\n")
	print.data.frame(round(x$dat,2))

}

