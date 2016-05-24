make.classes <- function(array, n) {
	array<-array[!is.na(array)]
	breaks<-seq(min(array),max(array), len=n)
	interval<-findInterval(array,breaks,rightmost.closed=T)
	counts<-tabulate(interval)
	counts<-counts[order(breaks[2:n], decreasing=T)]
	breaks<-breaks[order(breaks, decreasing=T)]
	counts<-counts/sum(counts)
	return(data.matrix(cbind(breaks,c(0,counts))))
}

