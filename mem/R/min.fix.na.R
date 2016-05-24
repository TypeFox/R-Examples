min.fix.na <-
function(minputdata){
	minputdata[minputdata==Inf]<-NA
	minputdata[minputdata==-Inf]<-NA
	if (sum(!is.na(minputdata))==0){
		return(NA)
	}else{
		return(min(minputdata,na.rm=T))
	}
}
