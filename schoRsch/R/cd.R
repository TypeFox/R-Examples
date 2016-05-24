cd <- function(x) {
	if (x=="..") {
		# Compute directory length
		dirs <- unlist(strsplit(getwd(),"/"))
		dirlgth <- nchar(dirs[length(dirs)])
		# Remove current directory from wd string
		newwd <- strtrim(getwd(),nchar(getwd())-dirlgth-1)
		a<-1
	} else {
		newwd <- paste(getwd(),"/",x,sep="")
	}
	setwd(newwd)
	getwd()
}