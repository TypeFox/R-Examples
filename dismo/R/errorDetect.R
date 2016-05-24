# Author: Jacob van Etten
# License GPL3
# Version 0.1
# October 2008



#Detect conversion error (degrees-minutes to decimals) 
.conversionError <- function(xy, min_ndigits){

	ndigits <- function(x){
		result <- pmax(0, nchar(abs(x))-nchar(abs(trunc(x, digits=0)))-1)
		return(result)
	}

	xy <- stats::na.omit(xy)
	number_digits <- cbind(ndigits(xy[,1]), ndigits(xy[,2]))
	index <- which((pmin(number_digits[,1], number_digits[,2]) + (as.numeric(abs(number_digits[,1]-number_digits[,2]))==1))>=min_ndigits)
	xy <- xy[index,]
	xy <- abs(xy)
	xy_dec <- xy - trunc(xy, digits=0)
	fr <- vector(length=10)
	for (i in 1:10){
	fr[i] <- length(subset(xy_dec, xy_dec[,1] > (i-1)/10 & xy_dec[,2]<i/10)[,1])}
	p.value <- chisq.test(fr, p=rep(0.1,times=10))$p.value
	return(list(p.value = p.value, frequencies = fr))
}
	
	