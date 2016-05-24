max.fix.na <-
function(minputdata){
	minputdata[minputdata==Inf]<-NA
	minputdata[minputdata==-Inf]<-NA
	if (sum(!is.na(minputdata))==0){
		return(NA)
	}else{
		return(max(minputdata,na.rm=T))
	}
}
