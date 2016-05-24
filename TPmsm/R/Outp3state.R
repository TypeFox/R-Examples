Outp3state <- function(data, v, ...) {
	return( data.frame("times1"=data[,v[1]], "delta"=as.integer(data[,v[3]] > data[,v[1]]), "times2"=data[,v[3]]-data[,v[1]], "time"=data[,v[3]], "status"=data[,v[4]], data[,-v]) )
}
