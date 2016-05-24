`read.production` <- function(data){
	read.zoo(data, sep=",", FUN=fmt.chron, header=TRUE)
	}