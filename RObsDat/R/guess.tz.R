guess.tz <- function(tz){
	stopifnot(!is.null(tz))
	return(switch(as.character(tz), "0"="GMT", stop(paste("undefined tz:", tz))))

}
