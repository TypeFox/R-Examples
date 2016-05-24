# Author: Jacob van Etten
# License GPL3
# Version 0.1
# October 2008


.precision <- function(xy){
	
	ndigits <- function(coord){
		ndig <- pmax(0,nchar(abs(coord))-nchar(abs(trunc(coord,digits=0)))-1)
		return(ndig)
	}
	
	seconds <- function(coord){
		decimals <- 0:60 / 60
		diffnearestfraction <- sapply(coord, function(x) min(abs(x - decimals)))
		precision <- 10^-ndigits(coord) <= diffnearestfraction
		return(precision)
	}
	
	le <- length(xy[,1])
	xy <- cbind(1:length(xy[,1]),xy)
	xy <- stats::na.omit(xy)
	index <- rep(1,times=length(xy[,1]))
	index[which(pmin(cbind(ndigits(xy[,2]),ndigits(xy[,3]))) > 0)] <- 2
	index[index == 2] <- index[index == 2] + as.numeric(seconds(xy[which(index == 2),2] * seconds(xy[which(index ==2),3])))
	xy <- stats::na.omit(xy)
	xy.dec <- xy - trunc(xy,digits=0)
	precision <- seconds(as.vector(xy.dec))
	precision <- pmax(cbind(precision[1:length(xy[,1])],precision[(length(xy[,1])+1):length(precision)]))
	index.final <- rep(0, times=le)
	index.final[xy[,1]] <- index
	return(index.final)
}