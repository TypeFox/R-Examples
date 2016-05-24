mwd.thrsh<-function(data,verbosity=3,...){
	if(try(require(wavethresh,quietly=TRUE))){
		pad=length(rep(0,2^ceiling(log2(length(data)))-length(data)))
		before=floor(pad/2)
		after=ceiling(pad/2)
		adj=c(rep(0,before),data,rep(0,after))
		data<-wavethresh::mwd(adj,verbose=ifelse(verbosity>1,
			TRUE,FALSE),...)
		data<-wavethresh::threshold(data,verbose=ifelse(verbosity>0,
			TRUE,FALSE),...)
		data<-wavethresh::wr(data,verbose=ifelse(verbosity>2,
			TRUE,FALSE),...)
		data<-data[-(1:before)]
		data<-data[-((length(data)-after+1):length(data))]
		return(data)
	}else{
		stop("package wavethresh not available.\n")
	}
}

